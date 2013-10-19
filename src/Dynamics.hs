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
               ,dynamicVerlet
               ,initializeThermo
               ,landauZener
               ,moveCoord
               ,moveVel
               ,noseHoover1
               ,noseHoover2
               )where

import qualified Data.List as DL
import qualified Data.Vector.Unboxed as VU
import Control.Applicative
import Control.Arrow ((&&&))
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
dynamicVerlet :: Molecule -> DT -> Energy -> Job -> String -> IO Molecule
dynamicVerlet !mol !dt totalEnergy job project = do
  let step1 = moveVel dt $ moveCoord dt mol
  step2 <- interactWith job project step1
  let newMol = scaleVel totalEnergy . moveVel dt $ step2
  landauZener mol newMol dt totalEnergy job project
  
{- |function to keep constant the enegy, scaling the kinetic energy to the difference
between the initial energy (the potential energy in the Frank-Condon Point) and
the current potential energy -}
scaleVel :: Energy -> Molecule -> Molecule
scaleVel totalEnergy mol = mol{getVel = scale (getVel mol) fac}
  where fac = if deltaE < 1.0e-5 then 1.0 else sqrt $ deltaE / ek
        deltaE = abs $ ( elecE mol) - totalEnergy
        (vs,ms) = getVel &&& getMass $ mol
        ek = calcEk vs ms                      
 
-- | Step to advance the coordinates
moveCoord :: DT -> Molecule -> Molecule
moveCoord dt mol@(Molecule xs vs fs ms _at _es _derv _st _fc _ws) = mol{getCoord = newCoord}
  where dt2   = dt * 0.5
        dtsq2 = dt * dt2
        fun sh = (! sh) `fmap` [xs, vs, fs]
        newCoord = computeUnboxedS $ fromFunction (extent xs) $
                   (\sh@(Z :.i ) ->  let [x,v,f] = fun sh
                                         m = ms ! (Z :. i `div` 3)
                                     in x + dt*v + dtsq2*f/m )

-- | Step to semi-advance the velocity, here the velocity is only integrated a half of DT
moveVel ::  DT -> Molecule -> Molecule
moveVel dt mol@(Molecule _xs vs fs  ms _at _es _derv _st _fc _ws) = mol{getVel = newVel}
  where dt2   = dt * 0.5
        fun sh = (! sh) `fmap` [vs, fs]
        newVel = computeUnboxedS $ fromFunction (extent vs) $
                   (\sh@(Z :.i ) ->
                      let [v,f] = fun sh
                          m = ms ! (Z :. i `div` 3)
                      in v + dt2*f/m )

        
-- =================> NOSÉ-HOOVER <==================

-- |The record thermostat is in charge of the bookkeeping of a chain of thermostat.
-- |In this case there are only 2 thermostats
bath :: Molecule -> DT -> Temperature -> Thermo -> (Molecule,Thermo)
bath mol dt t thermo@(Thermo q1 q2 vx1 vx2) =
  let (vel,mass) = getVel &&& getMass $ mol
      n = (\(Z:. m) -> fromIntegral m) $ extent vel
      ek = calcEk vel mass
      g2 = (q1*vx1^2 - t*kb)/q2
      g1  = (2.0*ek -3.0*n*t*CTES.kb)/q1
      vx4 = vx2 + g2*dt4
      vx3 = (\x-> x*exp(-vx4*dt8)) . (\x -> x + g1*dt4) . (\x -> x*exp(-vx4*dt8)) $ vx1

      s = exp(-vx3*dt2)
      newVel = R.computeS . R.map (*s) $ vel
      ek2 = calcEk newVel mass

      g3 = (2.0*ek2 -3.0*n*t)/q1
      vx5 = (\x -> x*exp (-vx4*dt8)) . (\x -> x + g3*dt4) $ vx3
      g4 = (q1*vx5 - t*kb) / q2
      vx6 = vx4 + g4* dt4

  in (mol{getVel = newVel}, thermo{thVx1 = vx5, thVx2 = vx6})

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

--  let  newThermo = initializeThermo (length . getAtoms $ mol) t 2.2e1



dynamicNoseHoover ::  Molecule -> DT -> Temperature -> Thermo -> Job -> String -> IO (Molecule,Thermo)  
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
appliedForce [m,n] modForceExt mol = mol {getForce = newForce}
  where [v1,v2] = fmap (\j -> let i = 3*(pred j) in fmap (forceVector !) [(Z:.i),(Z:.i+1),(Z:.i+2)]) [m,n]
        forceVector = getForce mol
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

-- =================> Hoppping Algortihm         
                      
-- | Naive hopping algorithm 
landauZener :: Molecule -> Molecule -> DT -> Energy-> Job -> String -> IO Molecule
landauZener old new dt totalEnergy job project = undefined

        

-- ===================> UTILITIES <==================

calcTotalEnergy :: Molecule -> Double
calcTotalEnergy mol = kinetic + potential
  where kinetic = calcEk vs ms
        (vs,ms) = getVel &&& getMass $ mol
        state = getElecSt mol
        currentEnergies = head . getEnergy $ mol
        potential = currentEnergies !! (calcElectSt mol)                         
                         
hopDown :: Molecule -> Molecule
-- hopDown mol = let st = getElecSt mol in mol{getElecSt = pred st}                      
hopDown = undefined

hopUp :: Molecule -> Molecule
--hopUp mol = let st = getElecSt mol in mol{getElecSt = succ st}
hopUp = undefined
                 

calcEk :: Array U DIM1 Double -> Array U DIM1 Double -> Double
calcEk vs ms = R.sumAllS . computeUnboxedS . fromFunction (extent vs) $
    (\sh@(Z:.i) -> let m = ms ! (Z :. i`div` 3)
                       v = vs ! sh
                   in  0.5*m*v^2)
                   
dotP ::  Array U DIM1 Double -> Array U DIM1 Double -> Double
dotP v1 v2 = sumAllS . computeUnboxedS $ R.zipWith (*) v1 v2

scale :: Array U DIM1 Double -> Double -> Array U DIM1 Double
scale v s = computeUnboxedS . R.map (*s) $ v    

elecE :: Molecule -> Energy
elecE = undefined 
-- elecE mol = es !! k
--   where (st,es) = getElecSt &&& getEner $ mol
--         k = fromEnum st

normalize :: [Double] -> [Double]
normalize xs = fmap (*(recip norm)) xs
  where norm = sqrt . sum $ DL.zipWith (*) xs xs
                   