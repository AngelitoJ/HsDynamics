{-# Language BangPatterns #-}


module InitialConditions {- (
                           initialConditions
                          ,initializeMolcasOntheFly
                          ,initializeMolcasTinker
                          ,initializeRC
                          ,initializeSystemOnTheFly
                          ,parseConnections
                          ,parseFileInput
                          ,parseInput
                          )-} where

import qualified Data.List as DL
import qualified Data.Vector.Unboxed as VU
import Control.Applicative hiding ((<|>))
import Control.Arrow ((&&&),first)
import Control.Lens 
import Control.Monad ((<=<),liftM,mplus)
import Data.Array.Repa as R
import Data.Array.Repa.Unsafe as R
import qualified Data.Vector as V
import Data.Char(toUpper,toLower)
import Data.Either (either)
import Data.Maybe (fromMaybe)
import qualified Data.Map as M
import qualified Data.Monoid as DM
import Data.Random.Normal
import System.IO.Strict as IOS
import System.Random
import Text.Parsec
import Text.Parsec.ByteString


-- internal modules
import APIparser
import CommonTypes
import Constants
import Gaussian
import Molcas
import ReactionCoordinate
import ParsecNumbers
import ParsecText
import TinkerQMMM

-- =====================> <====================



-- =====================> <===================

initializeRC :: Int -> M.Map Int (Int,Internals)
initializeRC numat = M.fromList [(i,(0,v0)) | i <- [0..39]]
  where dim = 3*numat-6
        v0 = VU.replicate dim 0.0

-- =========================> <========================

-- | Parser to read the initial conditions declare in the input file
parseFileInput :: MyParser () InitialDynamics -> FilePath ->  IO InitialDynamics 
parseFileInput parser fileName = do
  r <- parseFromFile parser fileName 
  case r of
       Left msg -> fail $ show msg
       Right xs -> return xs

parseInput :: MyParser () InitialDynamics
parseInput = do 
  elecSt          <- parserElectronicState
  time            <- parserTotalTime
  step            <- parserIntegrationStep
  (theory,basis)  <- parserTheory
  extForce        <- parserExtenalforce
  anchor          <- parserAnchor        -- atoms on which the external forces are applied
  project         <- parserProject  
  return $ InitialDynamics elecSt time step theory basis extForce anchor project


-- | initial Electronic state
parserElectronicState :: MyParser () Singlet
parserElectronicState = do
  keyword <- parseKeyword
  st <- (\x -> read x :: Singlet) <$> manyTill anyChar space
  anyLine
  if (keyword)  == "initialstate" 
         then return st
         else fail "An initial electronic state is mandatory"
                               
parserTotalTime :: MyParser () Double 
parserTotalTime = do
  keyword <- parseKeyword
  time    <- realNumber
  anyLine
  if keyword == "time" then return time
                       else fail " The total time of the dynamics is mandatory"
                       
parserIntegrationStep :: MyParser () Double
parserIntegrationStep = do
  keyword <- parseKeyword
  step    <- realNumber
  anyLine
  if keyword == "step" then return step
                       else fail " The total Integatin step of the dynamics is mandatory"
                       
-- | in case of Gaussian Job a level of theory is required                       
parserTheory :: MyParser () (TheoryLevel,String)
parserTheory = option (Unspecified,"Unspecified") $ do
  spaces 
  try (string "Theory") <|> try (string "theory")
  spaces >> char '=' >> spaces
  theory  <- parseLevel
  basis   <- manyTill anyChar space
  anyLine
  return (theory,basis)
  
  where parseLevel = try parserCasMolcas <|> try parserHF <?> "theory level"
        parserHF   = do
                    try (string "HF") <|>  try (string "hf")
                    return HF
        parserCasMolcas = do          
                   try (string "CASSCF") <|>  try (string "casscf") 
                   char '('
                   electrons <- intNumber
                   char ','
                   orbitals <- intNumber
                   manyTill anyChar (char '=')
                   rlxRoot <- intNumber       
                   rest <- manyTill anyChar (char ')')                   
                   return $ CASSCF (electrons,orbitals) rlxRoot  rest               
        funParser s p = case runP p () "" s of
                             Left msg -> error $ show msg
                             Right x  -> x  
  

parserExtenalforce :: MyParser () Double
parserExtenalforce  = option 0 $ try $ do
  spaces
  try (string "Force") <|> try (string "force")
  parseEqual
  force   <- realNumber
  anyLine
  return force
  
-- | a list of atom numbers where the external force is applied  
parserAnchor :: MyParser () [Int]
parserAnchor  = option [] $ try $ do
  spaces
  try (string "Anchor") <|> try (string "anchor")
  parseEqual
  anchor  <- (\xs -> read xs :: [Int]) <$> manyTill anyChar space
  anyLine
  return anchor

-- | A string which gives name to the job  
parserProject :: MyParser () String
parserProject = option [] $ try $ do
  spaces
  try (string "Project") <|> try (string "project")
  parseEqual
  project <- manyTill anyChar space
  anyLine
  return project
                       
parseKeyword :: MyParser () String
parseKeyword = do
  spaces
  keyword <- manyTill anyChar space
  parseEqual
  return $ toLower <$> keyword
  
parseEqual :: MyParser () ()
parseEqual = spaces *> char '=' *> spaces

-- =====================> Initial Conditions <====================
 
-- function to initialize on the fly molecular dynamics using Molcas  
initializeMolcasOntheFly :: FilePath -> Singlet -> Temperature -> IO Molecule
initializeMolcasOntheFly xyz st temp = do
  r <- parseFromFile parseMoleculeXYZ xyz
  case r of
       Left msg -> error $  show msg
       Right xs -> initializeBaseOnXYZ xs st temp
   
initializeMolcasZeroVel ::  FilePath -> Singlet -> Temperature -> IO Molecule
initializeMolcasZeroVel xyz st temp = do
  r <- parseFromFile parseMoleculeXYZ xyz
  case r of
       Left msg -> error $  show msg
       Right xs -> do 
                  mol <- initializeBaseOnXYZ xs st temp
                  let dim = 3 * (mol ^. getAtoms . to length)
                      zeroVel = R.fromUnboxed (Z:. dim) $ VU.replicate dim 0
                  return $ set getVel zeroVel mol    
                       
 
initializeBaseOnXYZ :: [(Label,[Double])] -> Singlet -> Temperature -> IO Molecule
initializeBaseOnXYZ xs st temp = do
  let (labels,coord) = unzip xs
      f2U (c:cs)    = toUpper c : cs
      numat         = DL.length labels
      dim           = 3*numat
      repaCoord     = fromListUnboxed (Z :. dim) . fmap (/a0) $ concat coord
      masses        = fmap (\k -> fromMaybe msg $ M.lookup (f2U k) atom2Mass) labels
      msg           = error "Sorry boy, but I don't know some of your input atoms"
      aumasses      = fromListUnboxed (Z:.numat) $ fmap (*amu) masses   
      forces        = R.fromUnboxed (Z:.dim) $ VU.replicate dim 0
  velocities <- genMaxwellBoltzmann aumasses temp 
  return $ defaultMol & getCoord .~ repaCoord 
                      & getVel   .~ velocities 
                      & getForce .~ forces
                      & getMass  .~ aumasses 
                      & getAtoms .~ labels 
                      & getElecSt.~ (Left st) 
                     
initializeMolcasTinker :: FilePath -> Singlet -> Temperature -> Int -> IO Molecule  
initializeMolcasTinker molcasInput st temp numat = do
  molcasQM   <- parserInputMolcasQM molcasInput $ parserGatewayQM numat
  let atoms  = takeLabelCoordinates <$> molcasQM
  initializeBaseOnXYZ atoms st temp
                                             
  where takeLabelCoordinates (AtomQM label xyz) = (label,xyz)  
  
-- |collects the initial required to initialize a dynamic on the fly
initializeSystemOnTheFly :: FilePath -> Singlet -> Temperature-> IO Molecule
initializeSystemOnTheFly file st temp = do
     let keywords = ["Coordinates","Grad","Masses","Charges"]
     pairs <- takeInfo keywords <=< parseGaussianCheckpoint $ file
     let getData = lookupLabel pairs
         [coord,grad,masses] = fmap ((\ys -> R.fromListUnboxed (Z:. DL.length ys) ys) . getData) ["Coordinates","Grad","Masses"]
         labels = fmap (charge2Label .floor) $ getData "Charges"
         aumasses = computeUnboxedS $ R.map (*amu) masses
         forces = computeUnboxedS $ R.map (negate) grad
     velocities <- genMaxwellBoltzmann aumasses temp
     return $ defaultMol & getCoord .~ coord 
                         & getVel   .~ velocities 
                         & getForce .~ forces
                         & getMass  .~ aumasses 
                         & getAtoms .~ labels 
                         & getElecSt.~ (Left st)                          
-- =============> <================
        
initialConditions :: [Double] -> [Double] -> [EnergyDerivatives] -> Connections -> Temperature -> IO Molecule
initialConditions !cartCoord !masses !dervEs !conex !t = do
     let numat = length masses
     velocities <- genMaxwellBoltzmann (R.fromListUnboxed (Z:. numat) masses) t
     let forces = R.fromUnboxed (Z:. numat*3) $ VU.generate (3*numat) (const 0)
         repaMass  = R.fromListUnboxed (Z:. numat) masses
         repaCoord = R.fromListUnboxed (Z:. numat*3) cartCoord
         internals = calculateInternals conex (R.toUnboxed repaCoord)
         initialMol = defaultMol
         fc = FC internals conex
     return $ initialMol & getCoord .~ repaCoord
                         & getVel   .~ velocities 
                         & getForce .~ forces
                         & getMass  .~ repaMass
                         & getElecSt.~ (Left S0)
                         & getDervEn.~ dervEs
                         & getFCStruc.~fc

genMaxwellBoltzmann :: Array U DIM1 Double -> Temperature -> IO (Array U DIM1 Double)
genMaxwellBoltzmann !ms !t = do
  let stds = fmap (\mi -> sqrt $ (t*kb/mi)) $ R.toList ms
  xs <- mapM (\std -> (take 3) `liftM` (normalsIO' (0,std)))  stds  
  let vs = concat xs
  return $ R.fromListUnboxed (Z:. length vs) vs

parseConnections :: FilePath -> IO Connections
parseConnections file = do
  r <- parseFromFile parseInternals file
  case r of
       Left msg -> error $ show msg
       Right conex -> return conex
  
-- ===================> <==================
-- | String to Nuclear Charge
atom2Mass :: M.Map String Double 
atom2Mass = M.fromList [("H",1.00782504),("C",12.0),("N",14.0),("O",16.0)]
  


