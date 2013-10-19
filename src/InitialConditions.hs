{-# Language BangPatterns #-}


module InitialConditions  (
                           initialConditions
                          ,initializeMolcasOntheFly
                          ,initializeRC
                          ,initializeSystemOnTheFly
                          ,parseConnections
                          ,parseInputDynamics
                          ) where

import qualified Data.List as DL
import qualified Data.Vector.Unboxed as VU
import Control.Applicative
import Control.Arrow ((&&&))
import Control.Monad ((<=<),liftM,mplus)
import Data.Array.Repa as R
import Data.Array.Repa.Unsafe as R
import Data.Char(toUpper)
import Data.Either (either)
import Data.Maybe (fromMaybe)
import qualified Data.Map as M
import qualified Data.Monoid as DM
import Data.Random.Normal
import System.IO.Strict as IOS
import System.Random
import Text.ParserCombinators.Parsec


-- internal modules
import APIparser
import CommonTypes
import Constants
import Gaussian
import ReactionCoordinate

-- =====================> <====================



-- =====================> <===================

initializeRC :: Int -> M.Map Int (Int,Internals)
initializeRC numat = M.fromList [(i,(0,v0)) | i <- [0..39]]
  where dim = 3*numat-6
        v0 = VU.replicate dim 0.0

-- =====> <==========

parseInputDynamics :: FilePath -> IO InitialDynamics
parseInputDynamics file = do
  r <- lines `fmap` IOS.readFile file
  let (ls,ts) = splitAt 5 r
      [st,t,dt,mf,a] = fmap (head .words) ls
      [time,delta,modForce] = fmap (\x -> read x :: Double)  [t,dt,mf]
      anchor = (read a) :: [Int]
      project = head . words $ head ts
  return $ InitialDynamics (read st) time delta modForce anchor project

-- function to initialize on the fly molecular dynamics using Molcas  
initializeMolcasOntheFly :: FilePath -> Singlet -> Temperature -> IO Molecule
initializeMolcasOntheFly xyz st temp = do
  dat <- IOS.readFile xyz
  case runParser parseMoleculeXYZ () "" dat of
       Left msg -> error $  show msg
       Right xs -> initializeBaseOnXYZ xs st temp
-- 
       
initializeBaseOnXYZ :: [(Label,[Double])] -> Singlet -> Temperature -> IO Molecule
initializeBaseOnXYZ xs st temp = do
  let (atoms,coord) = unzip xs
      f2U (c:cs)    = toUpper c : cs
      numat         = DL.length atoms
      dim           = 3*numat
      repaCoord     = fromListUnboxed (Z :. dim) . fmap (/a0) $ concat coord
      masses        = fmap (\k -> fromMaybe msg $ M.lookup (f2U k) atom2Mass) atoms
      msg           = error "Sorry boy, but I don't know some of your input atoms"
      aumasses      = fromListUnboxed (Z:.numat) $ fmap (*amu) masses     
  velocities <- genMaxwellBoltzmann aumasses temp 
  return $ defaultMol{
                     getCoord=repaCoord
                    ,getVel=velocities
                    ,getMass=aumasses
                    ,getAtoms=atoms
                    ,getElecSt=Left st
                     }
  
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
         initialMol = defaultMol
     velocities <- genMaxwellBoltzmann aumasses temp
     return $ initialMol{getCoord=coord,getVel=velocities,getForce=forces,
                         getMass=aumasses,getAtoms=labels,getElecSt=Left st}
                          

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
     return $ initialMol{getCoord=repaCoord,getVel=velocities,getForce=forces,
                         getMass=repaMass,getDervEn=dervEs,getElecSt = Left S0, getFCStruc= fc}


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
  
-- readArrayDIM1FIle :: FilePath -> IO (Array U DIM1 Double)
-- readArrayDIM1FIle file = do
--    vs <- readVectorFromFile file
--    return $ R.fromUnboxed (Z:. VU.length vs) vs
-- 
-- readArrayDIM2FIle :: FilePath -> IO (Array U DIM2 Double)
-- readArrayDIM2FIle file = do
--    vs <- readVectorFromFile file
--    let dim = VU.length vs
--    return $ R.fromUnboxed (Z:. dim :. dim) vs
--    
--   
-- readVectorFromFile :: FilePath -> IO (VU.Vector Double)
-- readVectorFromFile file = do  
--     s   <- L.readFile file
--     return $ parseL s
--  
-- -- Fill a new vector from a file containing a list of numbers.
-- parseL :: L.ByteString -> VU.Vector Double
-- parseL = VU.unfoldr step
--   where
--      step !s = case L.readDouble s of
--         Nothing       -> Nothing
--         Just (!k, !t) -> Just (k, L.tail t)  

