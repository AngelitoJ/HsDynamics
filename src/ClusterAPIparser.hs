{-# Language DeriveFunctor ,BangPatterns, TupleSections #-}

module APIparser where

import Control.Applicative
import Control.Arrow ((&&&))
import Control.Concurrent (threadDelay)
import Control.Concurrent.Async
import Control.Lens
import Control.Monad 
import Data.Array.Repa as R hiding ((++))
import Data.Complex
import Data.Either (lefts)
import Data.List as DL
import Data.Maybe (catMaybes,fromMaybe,isJust)
import qualified Data.Vector as V
import System.Directory 
import System.Exit(ExitCode( ExitFailure,ExitSuccess))
import System.Posix.User (getLoginName)
import System.Process (system,readProcess)
import Text.Parsec
import Text.Printf (printf)


-- Internal modules
import CommonTypes
import Constants
import Gaussian 
import Molcas 
import ParsecNumbers
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
     Molcas -> do
                 writeMolcasXYZ (project ++ ".xyz") mol
                 launchJob project
                 parseMolcas project ["Grad","Roots"] mol
                 
     MolcasTinker atomsQM tinkerPath -> do 
                          reWriteXYZtinker mol atomsQM tinkerPath project
                          launchTinker tinkerPath 
                          modifyTinkerNames project
                          launchMolcas project
                          parseMolcas project ["GradESPF"] mol
                                                                      
                 
     Gaussian theoryLevel -> do 
                            let input  = project ++ ".com"
                                format = "/apps/Chem/g09/formchk " ++ project ++ ".chk"
                                fchk   = project ++ ".fchk"        
                                out    = project ++ ".log"        
                            writeGaussJob theoryLevel project mol 
                            launchGaussian input
                            updateMultiStates out mol
                            
--      Palmeiro  -> do
       
                             
     Quadratic -> return $ calcgradQuadratic mol

-- Command to summit jobs in the Resmol cluster     
launchCluster :: Command -> Args -> IO ()
launchCluster cmd name = do
  [pid] <- lines `fmap` readProcess cmd [name] ""
  simpleLoop pid 10 not    -- has the calculation started?
  simpleLoop pid 20 id -- has the calculation finished?

-- | check the cluster queue each t seconds
simpleLoop :: String  -> Int -> (Bool -> Bool) -> IO ()
simpleLoop pid t fun = do
  jobsId <- (concatMap words . lines) <$>  readProcess "qstat" ["-r"] ""
  if (fun $ DL.elem pid jobsId)
        then threadDelay (10^6*t) >> simpleLoop pid t fun -- Haskell use miliseconds
        else return ()

-- | local jobs  
launchJob :: String -> IO ()
launchJob script = do
  exit <- async $ system $ script
  checkStatus <=< wait $ exit

-- | Job Status 
checkStatus :: ExitCode -> IO ()
checkStatus r =  case r of
       ExitSuccess ->  return ()
       ExitFailure _ -> putStrLn "Lord give us patience and resistence in the ass!!\n" >>
                        fail "Job Launch Failed (Bastard Fortranians!!)"
                        
-- ============> Shared Functions <=========
updateCoeffEnergies :: [Energy] -> [[Double]] -> Molecule -> Molecule
updateCoeffEnergies energies coeff mol = 
  let currentCoeff =  mol^.getCoeffCI
      currentEner  =  mol^.getEnergy
      (newCoeffs,newEnergy) = if length currentCoeff < 3 
         then (coeff : currentCoeff, energies : currentEner)
         else (coeff : init currentCoeff ,energies : init currentEner)
  in mol & getCoeffCI .~ newCoeffs 
         & getEnergy  .~ newEnergy
                                                
-- =======================> Tinker <===========                                                
launchTinker :: Command -> IO ()
launchTinker = undefined

-- | Tinker does not overwrite the .xyz, instead it writes thew new geometry optimization 
--   in a file ended in .xyz_2
modifyTinkerNames :: Project -> IO ()
modifyTinkerNames project = renameFile (root ++ "_2") root
  where root = project ++ ".xyz"
-- =======================> Molcas <============

launchMolcas :: Project ->  IO ()
launchMolcas project = launchCluster "Molcas" $ project ++ ".input"

parseMolcas :: Project -> [Label] -> Molecule -> IO Molecule
parseMolcas project labels mol = do
  pairs <- takeInfoMolcas labels <=< parseMolcasOutputFile $ project ++ ".out"
  let fun = lookupLabel pairs
      newMol = updateMolecule fun labels mol
  return newMol
                  
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
                                else  let [ModuleAuto _n _s modules] = filter Molcas.isAuto xs
                                      in Just modules     
                                           
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
getInfoFCHK :: [Label] -> FilePath -> Molecule -> IO Molecule      
getInfoFCHK keywords file mol = do
  pairs <- takeInfo keywords <=< parseGaussianCheckpoint $ file
  let grad = (\ys -> R.fromListUnboxed (Z:. DL.length ys) ys) $ lookupLabel pairs "Grad"
      forces = computeUnboxedS $ R.map (negate) grad
  return $ set getForce forces mol

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

parseInternals :: MyParser st (V.Vector InternalCoord)
parseInternals = do
  bonds     <- parseSection 2
  angles    <- parseSection 3
  dihedrals <- parseSection 4
  return $  bonds V.++ angles V.++ dihedrals
  
parseSection :: Int ->  MyParser st (V.Vector InternalCoord)
parseSection !n = do
  spaces
  nInternals <- intNumber
  skipLine
  conex <- count nInternals $ parseQ n
  return $ V.fromList conex
  
-- Since indexes start at zero the atom numbering should
-- be reduced by 1.

parseQ :: Int ->  MyParser st InternalCoord 
parseQ !n = do 
           xs <- count n (spaces >> intNumber)
           skipLine
           let getIndex x = pred (xs !! x)
           case n of 
                2 -> return $ Bond     (getIndex 0) (getIndex 1)
                3 -> return $ Angle    (getIndex 0) (getIndex 1) (getIndex 2)
                4 -> return $ Dihedral (getIndex 0) (getIndex 1) (getIndex 2) (getIndex 3)
     
     
-- =======================> ParseMolecules in XYZ format <=======================

parseMoleculeXYZ :: MyParser st [(Label,[Double])]  
parseMoleculeXYZ = do
    numat <- intNumber 
    count 2 skipLine
    geometry <- count numat parseAtoms
    return $ geometry 
    

parseAtoms :: MyParser st (Label,[Double])   
parseAtoms = do 
            spaces 
            label <- many1 alphaNum
            xs <- count 3 (spaces >> realNumber)
            skipLine
            return $ (label,xs) 
           

skipLine :: MyParser st ()
skipLine = manyTill anyChar newline >> return ()
     
            
-- =========> Utilities <=================
calcElectSt :: Molecule -> Int
calcElectSt mol = 
   case  mol^.getElecSt of
        Left s  -> fromEnum s
        Right s -> fromEnum s

        
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

writeMolcasXYZ :: FilePath -> Molecule -> IO ()
writeMolcasXYZ name mol = do
  let s = showAtoms mol
      numat = length $  mol^.getAtoms
      strAtoms  =  (show numat) ++ "\n"
      comment = "Angstrom\n"
      str = strAtoms ++ comment ++ s
  writeFile name str


writeGaussJob :: TheoryLevel -> String -> Molecule -> IO ()
writeGaussJob theoryLevel project mol =  do
  name <- getLoginName   
  let l1 = addNewLines 1 $ "%chk=" ++ project ++ ".chk"
      l2 = addNewLines 1 "%mem=4000Mb"
      l3 = addNewLines 1 "%nproc=2"
      l4 = addNewLines 1 $ "%scr=/scratch/" ++ name ++ "/"
      l5 = addNewLines 1 $ "%rwf=/scratch/" ++ name ++ "/"
      l6 = addNewLines 1 $ "#p " ++ theoryLevel ++ "  force  iop(1/33=1) nosymm"
      l7 = addNewLines 2 $ "# SCF=(MaxCycle=300,conver=7)"
      l8 = addNewLines 2 "save the old Farts, use Fortran."
      l9 = addNewLines 1 "0 1"
      atoms = addNewLines 1 $ showAtoms mol
      weights = if mol^.getElecSt == Left S0 then "" else addNewLines 5 $ " 0.5       0.5"
      result = foldl1 (++) [l1,l2,l3,l4,l5,l6,l7,l8,l9,atoms,weights]
  writeFile (project ++ ".com") result 
       
         
printMol :: Molecule -> String -> IO ()
printMol mol msg = appendFile "geometry.out" $ numat ++ (showAtoms mol)
  where numat = (show . length $ labels) ++ "\n" ++ msg ++ "\n"
        labels = mol^.getAtoms 

printData :: Molecule -> Int -> IO ()
printData mol step = do
  let st =  mol^.getElecSt 
      l1 = addNewLines 1 $ "step: " ++ (show step) 
      l2 = addNewLines 1 $ "electronic State: " ++ (show st)
      l3 = addNewLines 2 $ "potential energies: " ++ (concatMap (printf "%.6f  ") . head $  mol^.getEnergy)          
  appendFile "result.out" $ foldl1' (++) [l1,l2,l3]  
        
showAtoms :: Molecule -> String
showAtoms mol  = concatMap fun atoms
  where labels = mol^.getAtoms
        qs =  mol^.getCoord
        xs =  R.toList . computeUnboxedS . R.map (*a0) $  qs         
        atoms = zip labels [[xs !! (3*i), xs !! (3*i+1), xs !! (3*i+2)] | i <- [0..pred (length xs `div` 3)]]
        fun = \(s,[x,y,z]) -> printf "%s  %.5f  %.5f  %.5f\n" s x y z
         
addNewLines :: Int -> String -> String
addNewLines n s = let f = (++"\n")
                  in iterate f s !! n        
        
        
