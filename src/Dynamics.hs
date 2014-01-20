{-# Language FlexibleContexts,BangPatterns,ConstraintKinds #-}

module Dynamics (
                DT
               ,Energy
               ,Molecule(..)
               ,Thermo(..)
               ,Time
               ,calcEk
               ,calcElectSt
               ,calcTotalEnergy
               ,dynamicExternalForces
               ,dynamicNoseHoover
               ,initializeThermo
               ,moveCoord
               ,moveVel
               ,noseHoover1
               ,noseHoover2
               ,velocityVerletForces
               )where

import qualified Data.List as DL
import qualified Data.Vector.Unboxed as VU
import Control.Applicative
import Control.Arrow ((&&&))
import Control.Lens
import Control.Monad ((<=<),mplus)
import Data.Array.Repa as R
import Data.Array.Repa.Unsafe as R
import Data.Either (either)
import Data.Maybe (fromMaybe)

-- internal modules
import CommonTypes
import APIparser
import Constants as CTES 
               

-- ==========================> <==================

sumVector :: (VU.Unbox a, Num a) => Array U DIM1 a -> Array U DIM1 a -> Array U DIM1 a 
sumVector v1 v2 = computeUnboxedS $ v1 +^ v2

subVector :: (VU.Unbox a, Num a) => Array U DIM1 a -> Array U DIM1 a -> Array U DIM1 a 
subVector v1 v2 = computeUnboxedS $ v1 -^ v2

-- =========================> Velocity Verlet <==================================

{- |This is the constant energy function, Because we have an interphase with some
 external software in charge of the potential (like Molcas, Gaussian, etc)
 this function should be monadic-}
-- dynamicVerlet :: Molecule -> DT -> Energy -> Job -> String -> IO Molecule
-- dynamicVerlet !mol !dt totalEnergy job project = do
--   let step1 = moveVel dt $ moveCoord dt mol
--   step2 <- interactWith job project step1
--   scaleVel totalEnergy . moveVel dt $ step2
  
{- |function to keep constant the enegy, scaling the kinetic energy to the difference
between the initial energy (the potential energy in the Frank-Condon Point) and
the current potential energy -}
{-scaleVel :: Energy -> Molecule -> Molecule
scaleVel totalEnergy mol = over getVel (scale fac) mol
  where fac = if deltaE < 1.0e-5 then 1.0 else sqrt $ deltaE / ek
        deltaE = abs $ ( elecE mol) - totalEnergy
        ek = calcEk (mol^.getVel) (mol^.getMass)-}                      
 
-- | Step to advance the coordinates
moveCoord :: DT -> Molecule -> Molecule
moveCoord dt mol = set getCoord newCoord mol
  where [xs,vs,fs,ms] = fmap (mol ^.) [getCoord,getVel,getForce,getMass] 
        dt2   = dt * 0.5
        dtsq2 = dt * dt2
        fun sh = (! sh) `fmap` [xs, vs, fs]
        newCoord = computeUnboxedS $ fromFunction (extent xs) $
                   (\sh@(Z :.i ) ->  let [x,v,f] = fun sh
                                         m = ms ! (Z :. i `div` 3)
                                     in x + dt*v + dtsq2*f/m )

-- | Step to semi-advance the velocity, here the velocity is only integrated a half of DT
moveVel ::  DT -> Molecule -> Molecule
moveVel dt mol = set getVel newVel mol
  where [vs,fs,ms] = fmap (mol ^.) [getVel,getForce,getMass] 
        dt2   = dt * 0.5
        fun sh = (! sh) `fmap` [vs, fs]
        newVel = computeUnboxedS $ fromFunction (extent vs) $
                   (\sh@(Z :.i ) ->
                      let [v,f] = fun sh
                          m = ms ! (Z :. i `div` 3)
                      in v + dt2*f/m )

-- | Velocity-Verlet Algorithm 

velocityVerletForces ::  Molecule -> DT -> Job -> String -> Anchor -> Double -> IO Molecule  
velocityVerletForces !mol !dt job project anchor modForceExt = do  
  let  step1 = (\x y -> moveVel y . moveCoord y $ x) mol dt      
  step2 <- interactWith job project step1
  let newMol = appliedForce anchor modForceExt step2
  return $ moveVel dt newMol
          
-- =================> NOSÉ-HOOVER <==================

-- |The record thermostat is in charge of the bookkeeping of a chain of thermostat.
-- |In this case there are only 2 thermostats
bath :: Molecule -> DT -> Temperature -> Thermo -> (Molecule,Thermo)
bath mol dt t thermo@(Thermo q1 q2 vx1 vx2) =
  let vel  = mol ^. getVel
      mass = mol ^. getMass
      n    = mol ^. getAtoms . to (fromIntegral . length)
      ek   = calcEk vel mass
      g2   = (q1*vx1^2 - t*kb)/q2
      g1   = (2.0*ek -3.0*n*t*kb)/q1
      vx4  = vx2 + g2*dt4
      vx3  = (\x-> x*exp(-vx4*dt8)) . (\x -> x + g1*dt4) . (\x -> x*exp(-vx4*dt8)) $ vx1

      s = exp(-vx3*dt2)
      newVel = R.computeS . R.map (*s) $ vel
      ek2 = ek*s^2

      g3 = (2.0*ek2 -3.0*n*t*kb)/q1
      vx5 = (\x -> x*exp (-vx4*dt8)) . (\x -> x + g3*dt4) . (\x -> x*exp(-vx4*dt8)) $ vx3
      g4 = (q1*vx5^2 - t*kb) / q2
      vx6 = vx4 + g4* dt4

  in (set getVel newVel mol, Thermo q1 q2 vx5 vx6)

  where [dt2,dt4,dt8] = tail . take 4 . iterate (*0.5) $ dt
                 
  
-- | NoseHoover 1 and 2 advance the coordinates, while
-- | interchanging energy with the termostat
noseHoover1 :: Molecule ->  DT -> Temperature -> Thermo -> (Molecule,Thermo)
noseHoover1 mol dt t thermo = (newMol,newThermo)
  where (mol2,newThermo) = bath mol dt t thermo
        newMol = (\x y -> moveVel y . moveCoord y $ x) mol2 dt

noseHoover2 :: Molecule -> DT -> Temperature -> Thermo -> (Molecule,Thermo)
noseHoover2  mol dt t thermo = bath mol2 dt t thermo
  where mol2 = moveVel dt mol


dynamicNoseHoover :: Molecule -> DT -> Temperature -> Thermo -> Job -> String -> IO (Molecule,Thermo)  
dynamicNoseHoover !mol !dt !t thermo job project = do  
       let (step1,thermo1) = noseHoover1 mol dt t thermo      
       step2 <- interactWith job project step1
       return $ noseHoover2 step2 dt t thermo1
 
dynamicExternalForces ::  Molecule -> DT -> Temperature -> Thermo -> Job -> String -> Anchor -> Double -> IO (Molecule,Thermo)  
dynamicExternalForces !mol !dt !t thermo job project anchor modForceExt = do  
       let (step1,thermo1) = noseHoover1 mol dt t thermo      
       step2 <- interactWith job project step1
       let newMol = appliedForce anchor modForceExt step2
       return $ noseHoover2 newMol dt t thermo1
 
initializeThermo :: Int -> Temperature -> Thermo
initializeThermo numat t  = Thermo q1 q2 vx1 vx2
  where q1 = 3.0 * (fromIntegral numat)* t * kb / (freq^2)
        q2 = t * kb / freq^2
        vx1 = 0.0
        vx2 = 0.0
        freq = recip (2.2e1/au_time)


-- ================> EXTERNAL FORCES <======================

appliedForce :: Anchor -> Double -> Molecule -> Molecule
appliedForce xs modForceExt mol | null xs = mol 
appliedForce [m,n] modForceExt mol = set getForce newForce mol
  where [v1,v2] = fmap (\j -> let i = 3*(pred j) in fmap (forceVector !) [(Z:.i),(Z:.i+1),(Z:.i+2)]) [m,n]
        forceVector = mol ^. getForce
        modauN = modForceExt * (recip $ 10^9) / auN
        deltaF = fmap (*modauN) . normalize $ DL.zipWith (-) v2 v1
        forceExt = R.fromListUnboxed sh $ DL.concat [sparseList x m n deltaF | x <- [0..pred numat]]
        newForce = computeUnboxedS $ forceVector +^ forceExt
        numat = dim `div` 3
        sh@(Z:. dim) = extent forceVector

sparseList :: Int -> Int -> Int -> [Double] -> [Double]
sparseList i m n xs | i == (pred m) =xs
                    | i == (pred n) = fmap negate xs
                    | otherwise = take 3 . repeat $ 0.0

        
-- ===================> UTILITIES <==================

calcTotalEnergy :: Molecule -> Double
calcTotalEnergy mol = kinetic + potential
  where kinetic = calcEk (mol^.getVel) (mol^.getMass) 
        currentEnergies = mol ^. getEnergy . to head
        potential = currentEnergies !! (calcElectSt mol)                         
                         
                 

calcEk :: Array U DIM1 Double -> Array U DIM1 Double -> Double
calcEk vs ms = R.sumAllS . computeUnboxedS . fromFunction (extent vs) $
    (\sh@(Z:.i) -> let m = ms ! (Z :. i`div` 3)
                       v = vs ! sh
                   in  0.5*m*v^2)
                   
dotP ::  Array U DIM1 Double -> Array U DIM1 Double -> Double
dotP v1 v2 = sumAllS . computeUnboxedS $ R.zipWith (*) v1 v2

scale :: Double -> Array U DIM1 Double ->  Array U DIM1 Double
scale s v = computeUnboxedS . R.map (*s) $ v    

normalize :: [Double] -> [Double]
normalize xs = fmap (*(recip norm)) xs
  where norm = sqrt . sum $ DL.zipWith (*) xs xs
                   