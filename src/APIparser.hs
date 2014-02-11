{-# Language DeriveFunctor ,BangPatterns, TupleSections #-}

module APIparser where

import Control.Applicative
import Control.Arrow ((&&&),(***))
import Control.Concurrent (threadDelay)
import Control.Concurrent.Async
import Control.Lens
import Control.Monad 
import Data.Array.Repa as R hiding ((++))
import qualified Data.ByteString.Char8 as C
import qualified Data.ByteString.Lex.Double as L
import Data.Complex
import Data.Either (lefts)
import Data.List as DL
import Data.List.Split (chunksOf)
import Data.Maybe (catMaybes,fromMaybe,isJust)
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import System.Directory 
import System.Exit(ExitCode( ExitFailure,ExitSuccess))
import System.Posix.User (getLoginName)
import System.Process (system,readProcess)
import Text.Parsec
import Text.Parsec.ByteString (parseFromFile)
import Text.Printf (printf)


-- Internal modules
import CommonTypes
import Constants
import Gaussian 
import InternalCoordinates
import Logger
import Molcas 
import ParsecNumbers
import ParsecText
import QuadraticInterpolation
import TinkerQMMM
-- =======================> Data, Types and Instances <==========================

newtype ParseInfo a = ParseInfo {
                             runParseInfo :: [Maybe a] } deriving (Show,Functor)
                                                          
-- ==================> Instances <=========================


instance Monad ParseInfo where
  return a = ParseInfo $ [Just a]
  ma >>= f = let xs = runParseInfo ma
                 fun x = case x of
                              Nothing -> Nothing
                              Just w  -> g w                               
                 g = safeHeadMaybe . catMaybes . runParseInfo . f
                 ys = fmap fun xs
             in ParseInfo ys

safeHeadMaybe :: [a] -> Maybe a
safeHeadMaybe [] = Nothing
safeHeadMaybe (x:_) = Just x
-- 
instance MonadPlus ParseInfo where
  mzero = defaultInfo
  (ParseInfo xs) `mplus` (ParseInfo ys) = ParseInfo $ xs `mplus` ys
    
-- ===========> <============

defaultInfo :: ParseInfo a
defaultInfo = ParseInfo []

naturalTransf :: ParseInfo a -> [a]
naturalTransf (ParseInfo []) = []
naturalTransf (ParseInfo xs) = catMaybes xs


-- =================> Call External Programs <==================

interactWith :: Job -> String -> Molecule -> IO Molecule
interactWith job project mol = 
  case job of 
     Molcas inputData -> do
                 let io1 = writeMolcasXYZ (project ++ ".xyz") mol
                     io2 = writeFile (project ++ ".input") $ concatMap show inputData
                 concurrently io1 io2 
                 launchMolcasLocal project
                 parseMolcas project ["Grad","Roots"] mol
                 
     MolcasTinker inputData atomsQM molcasQM ->  do 
                 print "Update QM atoms"
                 reWriteXYZtinker mol atomsQM project
                 print "Launching tinker"
                 launchTinker project
                 print "rewrite Molcas input"
                 modifyMolcasInput inputData molcasQM project mol
                 print "launch Molcas"
                 launchMolcas project
                 print "Parsing Molcas output"
                 parseMolcas project ["GradESPF"] mol
                                                                      
                 
     Gaussian tupleTheory -> do 
                 let input  = project ++ ".com"
                     out    = project ++ ".log"        
                 writeGaussJob tupleTheory project mol 
                 launchJob $ "g09 " ++ input
                 print "Updating Roots info"
                 updateMultiStates out mol
                 
     GroundState tupleTheory -> do
                let [input,out,chk,fchk] = DL.zipWith (++) (repeat project) [".com",".log",".chk",".fchk"]
                writeGaussJob tupleTheory project mol 
                launchJob $ "g09 " ++ input
                launchJob $ "formchk " ++ chk
                getGradEnerFCHK fchk mol                 

                            
     Palmeiro conex dirs ->  launchPalmeiro conex dirs mol
       
                             
     Quadratic -> return $ calcgradQuadratic mol

-- Command to summit jobs in the Resmol cluster     
launchCluster :: Command -> Args -> IO ()
launchCluster cmd name = do
  [pid] <- lines `fmap` readProcess cmd [name] ""
  simpleLoop pid 10 not    -- has the calculation started?
  simpleLoop pid 20 id -- has the calculation finished?
  
  
--Command to summit jobs in the nodes of a cluster 
-- launchNode :: Command -> Args -> IO ()
-- launchNode cmd name = do
--  
-- | check the cluster queue each t seconds
simpleLoop :: String  -> Int -> (Bool -> Bool) -> IO ()
simpleLoop pid t fun = do
  jobsId <- (concatMap words . lines) <$>  readProcess "qstat" ["-r"] ""
  if (fun $ DL.elem pid jobsId)
        then threadDelay (10^6*t) >> simpleLoop pid t fun -- Haskell use miliseconds
        else return ()

-- | local jobs  
launchJob :: String -> IO ()
launchJob script = withAsync (system script) $ wait >=> checkStatus
  
-- | Job Status 
checkStatus :: ExitCode -> IO ()
checkStatus r =  case r of
       ExitSuccess   -> return ()
       ExitFailure _ -> putStrLn "Lord give us patience and resistence in the ass!!\n" >>
                        fail "Job Launch Failed (Bastard Fortranians!!)"
                        
parallelLaunch :: [IO a] -> IO [a]
parallelLaunch actions = foldr conc (return []) actions
  where conc io ios = do
             (x,xs) <- concurrently io ios
             return (x:xs)

                        
-- ============> Tully Updates <=========
updateCoeffEnergies :: [Energy] -> [[Double]] -> Molecule -> Molecule
updateCoeffEnergies energies coeff mol = 
  let currentCoeff =  mol^.getCoeffCI
      currentEner  =  mol^.getEnergy
      (newCoeffs,newEnergy) = if length currentCoeff < 3 
         then (coeff : currentCoeff, energies : currentEner)
         else (coeff : init currentCoeff ,energies : init currentEner)
  in mol & getCoeffCI .~ newCoeffs 
         & getEnergy  .~ newEnergy
                              
updateNewJobInput :: Job -> Molecule -> Job
updateNewJobInput job mol = case job of
                                 Molcas   theory -> updateMolcasInput theory rlxroot
                                 otherwise       -> job               
                                                              
  where rlxroot = succ $ calcElectSt mol                                                              

updateMolcasInput ::  [MolcasInput String] -> Int -> Job
updateMolcasInput xs rlxroot = Molcas $ fmap modifyInputRas xs
 where modifyInputRas x = case x of
                               (RasSCF _ h t) -> RasSCF rlxroot h t
                               otherwise      -> x
                
-- =======================> Tinker <===========                                                
launchTinker :: String -> IO ()
launchTinker project = do
   let root   ="/home/marco/7.8.dev/tinker-5.1.09/source/dynaqmmm.x  "
       name   = project ++ ".xyz"                    -- first argument : .xyz file
       suffix = "  10  0.1 0.01 298 0.734 0.723 > /dev/null 2>&1"  -- argument : -> step (1) -> dt (default 1) -> dump (default 0.1) -> temperature (298 K)
                                                                   -- -> scaling factor first HLA  -> scaling factor second HLA
       job = root ++ name ++ suffix
   launchJob job
        
        
-- ======================> Palmeiro <==============
launchPalmeiro :: Connections -> [FilePath] -> Molecule -> IO Molecule
launchPalmeiro conex dirs mol =  do
   let internals = calcInternals conex mol
       carts     =  mol ^. getCoord
   (e1,f1) <- interpolation conex internals carts (dirs !! 0)
   (e2,f2) <- interpolation conex internals carts (dirs !! 1)
   return $ mol & getForce  .~ f1 
                & getEnergy .~ [[e1,e2]]
  
-- Because Palmeiro set of Utilities required a Directory for each electronic state then
-- there are created as many directories as electronic states involved 
interpolation :: Connections -> Internals -> Coordinates -> FilePath -> IO (Energy,Force)  
interpolation conex qs carts dir = do
  pwd <- getCurrentDirectory
  setCurrentDirectory $ pwd ++ dir
  print $ "Current Directory" ++ pwd ++ dir
  writePalmeiroScript qs
  launchJob "gfortran  Tools.f90 Diag.f90 CartToInt2.f90 FitSurfMN_linux3.f90 PalmeiroScript.f90 -o PalmeiroScript.o -L/usr/lib64/ -llapack -lblas"
  launchJob "chmod u+x PalmeiroScript.o"
  launchJob "./PalmeiroScript.o > /dev/null 2>&1 "
  gradInter <- readVectorUnboxed "Gradient.out"
  energy    <- (head . VU.toList) `fmap` readVectorUnboxed "Energy.out"
  grad      <- transform2Cart conex gradInter carts 
  let force = R.computeUnboxedS $  R.map negate grad
  setCurrentDirectory pwd
  print gradInter
  return (energy,force)

parserFileInternasCtl :: FilePath -> IO Connections
parserFileInternasCtl name = do
  r <- parseFromFile parserCtl name
  case r of
       Left err  -> error $ show err
       Right xs  -> return xs
  
parserCtl ::  MyParser st (Connections)
parserCtl = do 
  manyTill anyChar $ try $ string "[InternalCoord]"
  anyLine
  parseInternals
   
processMolcasOutputFile fname = do
        putStrLn $ "Processing file:" ++ (show fname) ++ "\n"
        input  <- C.readFile fname
        case (runParser molcasOutParser defaultParserState fname input) of
             Left err  -> print err
             Right xs  -> mapM_ print xs   
   
getSuffixFile ::  FilePath -> String -> IO FilePath
getSuffixFile path suff = do
  xs <- (filter (`notElem` [".",".."])) <$> (getDirectoryContents path)
  let ctl = head $ filter (isSuffixOf suff) xs
  return $ if null ctl then error "no .ctl file found" 
		       else ctl
  
-- =======================> Molcas <============

launchMolcas :: Project ->  IO ()
launchMolcas project = launchCluster "Molcas" $ project ++ ".input"

launchMolcasLocal :: Project ->  IO ()
launchMolcasLocal project = do
   let input = project ++ ".input"
       out   = project ++ ".out"
       err   = project ++ ".err"
   launchJob $ "molcas  " ++ input ++ ">  " ++ out ++ " 2>  " ++ err 


modifyMolcasInput :: [MolcasInput String] -> [AtomQM] -> Project -> Molecule -> IO ()
modifyMolcasInput inputData molcasQM project mol = do
       let newInput = fmap fun inputData
           fun dat = case dat of
                       Gateway x -> writeGateway molcasQM mol
                       otherwise -> dat
       writeFile (project ++ ".input") $ concatMap show newInput
       
       
writeGateway :: [AtomQM] -> Molecule -> MolcasInput String
writeGateway atoms mol = Gateway $ concat $ DL.zipWith3 fun xs atoms [1..]

  where xs        = fmap printxyz $ chunksOf 3 qs
        printxyz [x,y,z] = printf "%12.5f  %.5f  %.5f" x y z
        qs        = mol ^. getCoord . to R.toList
        symbols   = mol ^. getAtoms
        between x = " Basis set\n" ++ x ++ " End of Basis\n"
        spaces    = (" "++)
        ans       = spaces . (++"     Angstrom\n")
        mm        = spaces . (++"...... / MM\n")
        fun :: String -> AtomQM -> Int -> String
        fun x (AtomQM label _xyz typo) i = 
                let num = label ++ (show i)
                in case typo of
                     QM basis -> between $  (spaces label)  ++ "." ++ basis ++ "\n" ++ (ans $ num ++ x)
                     MM -> between $ (mm label) ++ (ans $ num ++ x) ++ " Charge=  -0.000000\n"                                          
        
        
parseMolcas :: Project -> [Label] -> Molecule -> IO Molecule
parseMolcas project labels mol = do
  pairs <- takeInfoMolcas labels <=< parseMolcasOutputFile $ project ++ ".out"
  let fun = lookupLabel pairs
      newMol = updateMolecule fun labels mol
  return newMol

parseMolcasInput :: FilePath -> IO [MolcasInput String]
parseMolcasInput fname = do
  input  <- C.readFile fname
  case runParser parserMolcasInput defaultParserState fname input of
       Left msg -> error $ show msg
       Right xs -> return xs

  
updateMolecule :: (Label ->[Double]) -> [Label] -> Molecule -> Molecule
updateMolecule fun labels mol = foldl' step mol labels
  where step acc l = case l of
                         "Grad"     -> let forces = (\ys -> R.fromListUnboxed (Z:. length ys) $ fmap negate ys) $ fun "Grad"
                                       in set getForce forces acc
                         "Roots"    -> updateCIroots fun acc
                         "GradESPF" -> let forces = (\ys -> R.fromListUnboxed (Z:. length ys) $ fmap negate ys) $ fun "GradESPF"
                                       in set getForce forces acc
                         otherwise -> error "unknown Label"                    
                   
updateCIroots :: (Label -> [Double]) -> Molecule -> Molecule
updateCIroots f mol = 
   let energies = f "Energies"
       nroots   = length energies
       cis = [f $ "CIroot" ++ (show i) | i <- [0..pred nroots]]
   in updateCoeffEnergies energies cis mol
                                   
takeInfoMolcas :: [Label] -> Either ParseError [MolBlock] -> IO [(Label,[Double])]
takeInfoMolcas labels eitherInfo =
  case eitherInfo of
       Left s -> error . show $ s
       Right xs -> return . naturalTransf $ parseDataMolcas xs labels
       
parseDataMolcas :: [MolBlock] ->  [Label] -> ParseInfo (Label,[Double])
parseDataMolcas molcasBlocks keywords =  foldl1' (mplus) $ fmap lookLabel  keywords
  where moduleAuto = fromMaybe (error "Molcas calculation has failed") $ getModuleAuto molcasBlocks
        lookLabel key = foldl' parseMB defaultInfo moduleAuto
          where parseMB !acc !mb = let x = molcasParsers key mb
                                   in  acc `mplus` x
                   

molcasParsers :: Label -> ModuleData -> ParseInfo (Label,[Double]) 
molcasParsers key (ModuleData _name _string xs)  = 
  case key of 
       "Grad"     -> getAlaskaGrad xs 
       "Roots"    -> getRASCFCI xs
       "GradESPF" -> getGradESPF xs
       otherwise  -> error "Unkwown label"
                                         
getModuleAuto :: [MolBlock] ->  Maybe [ModuleData]
getModuleAuto xs  =  if null xs then Nothing
                                else  let x  = filter Molcas.isAuto xs
                                      in case x of 
                                              [ModuleAuto _n _s modules] -> Just modules
                                              otherwise -> Nothing
                                           
getAlaskaGrad :: [SectionData] -> ParseInfo (Label,[Double]) 
getAlaskaGrad xs = if null xs then ParseInfo [Nothing]
                              else let ys = filter Molcas.isMolGrads xs
                                   in case ys of
                                           [] -> ParseInfo [Nothing]
                                           [AlaskaMolecularGradients _label t _gradESPF] -> return $ ("Grad",tuples2List t)

getGradESPF :: [SectionData] -> ParseInfo (Label,[Double]) 
getGradESPF xs = if null xs then ParseInfo [Nothing]
                              else let ys = filter Molcas.isMolGrads xs
                                   in case ys of
                                           [] -> ParseInfo [Nothing]
                                           [AlaskaMolecularGradients _label _t gradESPF] -> return $ ("GradESPF", concat gradESPF)

 
getRASCFCI :: [SectionData] -> ParseInfo (Label,[Double])
getRASCFCI xs = if null xs then ParseInfo [Nothing]
                           else let ys = filter Molcas.isCICoeffs xs
                                in case ys of
                                        [] -> ParseInfo [Nothing]   
                                        otherwise ->  getRoots ys

getRoots :: [SectionData] -> ParseInfo (Label,[Double])
getRoots xs = let (es,cis) = unzip $ fmap fun xs
                  energies = return ("Energies",es)
                  nroots   = length es
                  coefficients = foldl1' (mplus) [return ("CIroot"++(show i),cis !! i) | i <- [0..pred nroots]]
              in energies `mplus` coefficients
                
  where fun (RasSCFCI _r e ci) = 
            let cfs = fmap (\(_number,_label,coeff,_weight) -> coeff) ci 
            in (e,cfs)
                                              
tuples2List :: [(String,Char,Double)] -> [Double]
tuples2List xs = fmap go xs
  where go (s,c,v) = v
        dim = DL.length xs 

-- ===========> Gaussian Interface <=========
-- |take info from the formated check point
getGradEnerFCHK :: FilePath -> Molecule -> IO Molecule      
getGradEnerFCHK file mol = do
  let keywords = ["Grad","Energy"]
  pairs <- takeInfo keywords <=< parseGaussianCheckpoint $ file
  let grad   = (\ys -> R.fromListUnboxed (Z:. DL.length ys) ys) $ lookupLabel pairs "Grad"
      energy = lookupLabel pairs "Energy"
      forces = computeUnboxedS $ R.map (negate) grad
  return $ mol & getForce  .~ forces 
               & getEnergy .~ [energy]
    
-- | Monadic Lookup function base on keywords  
takeInfo :: [Label] -> Either ParseError [GauBlock] -> IO [(Label,[Double])]
takeInfo labels eitherInfo =
  case eitherInfo of
       Left s   -> error . show $ s
       Right xs -> return . naturalTransf $ parseDataGaussian xs labels

                                     
updateMultiStates :: FilePath -> Molecule -> IO Molecule
updateMultiStates file mol = do
  r <- parseLogGaussian numat file
  case r of
       Left msg -> error . show $ msg
       Right (GaussLog xs gs) -> return $ updateGrads gs . updateCASSCF xs $ mol
       
  where numat = mol ^. getAtoms . to length
        sh = mol ^. getForce . to extent
        negateRepa = computeUnboxedS . R.map (negate) 
        updateGrads gs mol = let grads = fmap (negateRepa . fromListUnboxed sh) gs
                                 st = calcElectSt mol
                                 g = grads !! st
                             in set getForce g mol    
                                                              
updateCASSCF :: [EigenBLock] -> Molecule ->  Molecule
updateCASSCF xs mol =
      let energies = fmap getRootEnergy xs
          coeff = fmap getRootCoeff xs
      in updateCoeffEnergies energies coeff mol                                            
         
-- summit a job in the cluster queue                                      
launchGaussian :: Project -> IO ()
launchGaussian project = launchCluster "Launch_gaussian" project
 

parseDataGaussian :: [GauBlock] ->  [Label] -> ParseInfo (Label,[Double])
parseDataGaussian gaussBlocks keywords =  ParseInfo $ lookLabel `fmap` keywords
  where lookLabel key = (key,) `fmap` DL.foldl' parseGB Nothing gaussBlocks
            where parseGB !acc !gb =
                        case gb of
                             RGauBlock label n xs ->
                                        let x = gaussParsers key label xs
                                        in  acc `mplus` x
                             otherwise ->  acc                                                             

gaussParsers :: Label -> Label -> [Double] -> Maybe [Double]
gaussParsers x =
                case x of
                     "Coordinates" -> gaussCoord
                     "Energy"      -> gaussEnergy
                     "Grad"        -> gaussGrad
                     "Hess"        -> gaussHess
                     "Masses"      -> gaussMass
                     "Charges"     -> gaussElems
                     "Coeff"       -> gaussCoeff


gaussCoord :: String -> [Double] -> Maybe [Double]
gaussCoord s xs = if (words s) == ["Current","cartesian","coordinates"] then (Just xs) else Nothing

gaussEnergy :: String -> [Double] -> Maybe [Double]
gaussEnergy s xs = if (words s) == ["Total","Energy"] then (Just xs) else Nothing

gaussGrad  :: String -> [Double] -> Maybe [Double]
gaussGrad s xs = if (words s) == ["Cartesian","Gradient"] then (Just xs) else Nothing

gaussHess  :: String -> [Double] -> Maybe [Double]
gaussHess s xs = if (words s) == ["Cartesian","Force","Constants"] then (Just xs) else Nothing

gaussMass  :: String -> [Double] -> Maybe [Double]
gaussMass s xs = if (words s) == ["Real","atomic","weights"] then (Just xs) else Nothing

gaussElems  :: String -> [Double]-> Maybe [Double]
gaussElems s xs = if (words s) == ["Nuclear","charges"] then (Just xs) else Nothing

gaussCoeff :: String -> [Double]-> Maybe [Double]
gaussCoeff = undefined
                             
-- ================> Parser Internal Coordinates <===============

parserFileInternas :: FilePath -> IO Connections
parserFileInternas name = do
  r <- parseFromFile parseInternals name
  case r of
       Left err  -> error $ show err
       Right xs  -> return xs

parseInternals :: MyParser st (Connections)
parseInternals = do
  bonds     <- parseSection 2
  angles    <- parseSection 3
  dihedrals <- parseSection 4
  return $ V.concat [bonds, angles, dihedrals]
  
parseSection :: Int ->  MyParser st Connections 
parseSection !n = do
  spaces
  nInternals <- intNumber
  anyLine
  conex <- count nInternals $ parseQ n
  return $ V.fromList conex
  
-- Since indexes start at zero the atom numbering should
-- be reduced by 1.

parseQ :: Int ->  MyParser st InternalCoord 
parseQ !n = do 
           xs <- count n (spaces >> intNumber)
           anyLine
           let getIndex x = pred (xs !! x)  -- atomic index begins at 1  
           case n of 
                2 -> return $ Bond     (getIndex 0) (getIndex 1)
                3 -> return $ Angle    (getIndex 0) (getIndex 1) (getIndex 2)
                4 -> return $ Dihedral (getIndex 0) (getIndex 1) (getIndex 2) (getIndex 3)
     
     
-- =======================> ParseMolecules in XYZ format <=======================

parserGeomVel :: FilePath -> Molecule -> IO [Molecule]
parserGeomVel xyz mol = do
  r <- parseFromFile (many1 $ parseMol mol parseAtomsVel ) xyz
  case r of
       Left msg -> error $  show msg
       Right xs -> return xs


parseMoleculeXYZ :: MyParser st (Label,[Double]) -> MyParser st [(Label,[Double])]  
parseMoleculeXYZ parser = do
    numat <- intNumber 
    count 2 anyLine
    geometry <- count numat parser
    return $ geometry 

parseMol :: Molecule -> MyParser st (Label,[Double])  -> MyParser st Molecule
parseMol mol parser = do
    numat <- intNumber 
    anyLine
    es <- count 2 (spaces >> realNumber)
    anyLine 
    geometry <- count numat parser
    let new = updatePosMom geometry mol 
    return $ new & getEnergy .~ [es] 
    
    
parseAtoms :: MyParser st (Label,[Double])   
parseAtoms = do 
            spaces 
            label <- many1 alphaNum
            xs <- count 3 (spaces >> realNumber)
            anyLine
            return $ (label,xs) 
           
parseAtomsVel :: MyParser st (Label,[Double])   
parseAtomsVel = do 
            spaces 
            label <- many1 alphaNum
            xs    <- count 6 (spaces >> realNumber)
            anyLine
            return $ (label,xs)      
            
-- =========> Utilities <=================
calcElectSt :: Molecule -> Int
calcElectSt mol = 
   case  mol^.getElecSt of
        Left s  -> fromEnum s
        Right s -> fromEnum s


readVectorUnboxed :: FilePath -> IO (VU.Vector Double)
readVectorUnboxed file = C.readFile file >>= \s -> return $ parseUnboxed s

parseUnboxed ::C.ByteString -> VU.Vector Double
parseUnboxed = VU.unfoldr step
  where
    isspace x = if x == ' ' then True else False     
    step !s = let b = C.dropWhile isspace s
              in case L.readDouble b of
                 Nothing       -> Nothing
                 Just (!k, !t) -> Just (k, C.tail t)         
        
updatePosMom :: [(Label,[Double])] -> Molecule -> Molecule 
updatePosMom as mol = mol & getCoord .~ repaCoord 
                          & getVel   .~ velocities  
                          
  where repaCoord   = R.fromListUnboxed (Z :. dim) xs                         
        velocities = R.fromListUnboxed (Z :. dim) vs 
        dim        = 3 * length as
        (xs,vs)    = concat *** concat $ unzip . fmap (splitAt 3 . snd ) $ as
         
-- =============================> <===============================
printGnuplot :: MatrixCmplx -> Molecule -> IO ()
printGnuplot matrix mol = do
  let energies = mol ^. getEnergy . to head
      es = concatMap (printf "%.5f  ") energies
      st = printf " %.5f " (energies !! calcElectSt mol)
      p1 = printf " %.5f " . realPart . (^2) $ matrix !! 0
      p2 = printf " %.5f " . realPart . (^2) $ matrix !! 3
      s  = DL.foldl1' (++) ["gnuplot: ",es,st,p1,p2]
  appendFile "TullyOutput" $ s ++ "\n\n"        

        
-- ==============================> Printing Molecule and scripts <==============

writePalmeiroScript :: Internals -> IO ()
writePalmeiroScript qs = do
  let spaces = concat $ replicate 5 " "
      inpQ = spaces ++ (writeInternasFortran qs)
      l1 = concatMap (spaces++) ["Program FittingToSurface\n","Use MultiNodalFitSurf\n","Implicit None\n","Real(kind=8),allocatable :: Inpq(:),g(:),H(:)\n","Real(kind=8) :: E\n","Integer i,IErr,NumberCoord\n","Logical :: NewSimplex\n\n"]
      dim = (show $ VU.length qs)
      l2 = spaces ++ "NumberCoord="++ dim ++ "\n\n"
      l3 = concatMap (spaces++) ["Call ReadDataFiles('./',.True.,IErr)\n\n","Allocate(Inpq(1:NumberCoord),g(1:NumberCoord),H(1:NumberCoord*(NumberCoord+1)/2))\n"]      
      l4 = concatMap (spaces++) ["Call InterpolatePES(Inpq,IErr,NewSimplex,E,g,H)\n", "If (IErr>0) Print *,'InterpolatePES Has FAILED!.'\n", "Print *, 'Gradient'\n",
                                 "open (unit=112, file='Gradient.out',status='replace',action='write')\n",
                                 "open (unit=113, file='Energy.out',status='replace',action='write')\n",
                                 "write(unit=112,fmt=*)  (real(g(i)) ,i=1," ++ dim ++ ")\n",                                 
                                 "write(unit=113,fmt=*) (real(E))\n",
                                 "close(unit=112)\n",
                                 "close(unit=113)\n\n","End Program\n"]
  writeFile "PalmeiroScript.f90" $  l1 ++ l2 ++ l3 ++ inpQ ++ l4
      
writeInternasFortran :: Internals -> String
writeInternasFortran qs = 
  let {-xss = chunksOf 8 . init . concatMap (printf "%2.6f,") $ VU.toList qs-}
      spaces    = concat $ replicate 5 " "
      (hs:tss)  = chunksOf 8 . fmap (printf "%2.6f,") $ VU.toList qs
      blocks    = foldl' fun (DL.concat hs) tss
      fun acc v = acc ++ "&\n" ++ spaces ++ "&" ++ DL.concat v
  in "Inpq=(/" ++ (init blocks) ++ "/)\n\n"  -- init remove the last ","
  
  
writeMolcasXYZ :: FilePath -> Molecule -> IO ()
writeMolcasXYZ name mol = do
  let s = showCoord mol
      numat = length $  mol^.getAtoms
      strAtoms  =  (show numat) ++ "\n"
      comment = "Angstrom\n"
      str = strAtoms ++ comment ++ s
  writeFile name str


writeGaussJob :: (TheoryLevel,Basis)  -> String -> Molecule -> IO ()
writeGaussJob (theory,basis) project mol =  do
--   name <- getLoginName   
  let l1 = addNewLines 1 $ "%chk=" ++ project 
      l2 = addNewLines 1 "%mem=4000Mb"
      l3 = addNewLines 1 "%nproc=4"
--       l4 = addNewLines 1 $ "%scr=/scratch/" ++ name ++ "/"
--       l5 = addNewLines 1 $ "%rwf=/scratch/" ++ name ++ "/"
      l6 = addNewLines 1 $ "#p " ++ (show theory) ++ basis ++ "  force  iop(1/33=1) nosymm"
      l7 = addNewLines 2 $ "# SCF=(MaxCycle=300,conver=7)"
      l8 = addNewLines 2 "save the old Farts, use Fortran."
      l9 = addNewLines 1 "0 1"
      atoms = addNewLines 1 $ showCoord mol
      st    = calcElectSt mol 
      weights = addNewLines 5 $  if st == 0 then "" else " 0.5       0.5"
      result = foldl1 (++) [l1,l2,l3,l6,l7,l8,l9,atoms,weights]
  writeFile (project ++ ".com") result 
                           

printMol :: Molecule -> String -> Logger -> IO ()         
printMol mol msg logger = logMessage logger $ numat ++ (showPositionVel mol)
  where numat = (show . length $ labels) ++ "\n" ++ msg ++ "\n"
        labels = mol^.getAtoms 
          
printData :: Molecule -> Int -> Logger -> IO ()        
printData mol step logger = do
  let st =  mol^.getElecSt 
      l1 = addNewLines 1 $ "step: " ++ (show step) 
      l2 = addNewLines 1 $ "electronic State: " ++ (show st)
      l3 = addNewLines 2 $ "potential energies: " ++ (concatMap (printf "%.6f  ") . head $  mol^.getEnergy)          
  logMessage logger $ foldl1' (++) [l1,l2,l3]  
        
        
showCoord :: Molecule -> String
showCoord mol  = concat $ DL.zipWith fun labels xs 
  where labels = if null ls then repeat "" else ls
        ls =  mol^.getAtoms
        qs = mol^.getCoord
        xs = chunksOf 3 . R.toList . computeUnboxedS . R.map (*a0) $  qs
        fun l x  = (printf "%s" l) ++ (printxyz x) ++ "\n"
        printxyz [x,y,z] = printf "%12.5f  %.5f  %.5f" x y z
         
showPositionVel :: Molecule -> String
showPositionVel mol  = concat $ DL.zipWith3 fun labels xs vs
  where labels = mol^.getAtoms
        qs = mol^.getCoord
        xs = chunksOf 3 . R.toList . computeUnboxedS . R.map (*a0) $  qs
        vs = mol^.getVel . to (chunksOf 3 . R.toList)
        fun l x v = (printf "%s" l) ++ (printxyz x) ++ (printxyz v) ++ "\n"
        printxyz [x,y,z] = printf "%12.5f  %.5f  %.5f" x y z
                   
addNewLines :: Int -> String -> String
addNewLines n s = let f = (++"\n")
                  in iterate f s !! n        
        
        
