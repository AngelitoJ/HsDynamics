
module QuadraticInterpolation where

import Data.Array.Repa as R
import Data.Array.Repa.Unsafe as R
import qualified Data.Vector.Unboxed as VU


-- Internal Modules 
import CommonTypes
import LinearAlgebra 
import ReactionCoordinate


-- ==================> <==================
type DeltaQ = Array U DIM2 Double

type GradMtx = Array U DIM2 Double

type HessMtx = Array U DIM2 Double
-- =============> <=================
        
calcgradQuadratic :: Molecule -> Molecule
calcgradQuadratic mol =  undefined
{- let [state1,state2] = getDervEn mol
     EnergyDerivatives e1 grad1 hess1 = state1
     EnergyDerivatives e2 grad2 hess2 = state2
     (Z:.dim) = extent grad1
     [gmtx1,gmtx2] = fmap (\x -> R.computeUnboxedS $ reshape (Z:.dim :.1) x) [grad1,grad2]
     currentX = R.toUnboxed $ getCoord mol
     FC geomFC conex = getFCStruc mol
     currentQ = calculateInternalsFC geomFC conex currentX
     dq = R.fromUnboxed (Z:.dim :.1) $ VU.zipWith (-) currentQ geomFC     
     (energy1,gradQ1) = calcQuadraticEnergy dq e1 gmtx1 hess1 
     (energy2,_gradQ2) = calcQuadraticEnergy dq e2 gmtx2 hess2 
     newForce = VU.fromUnboxed (extent gradQ1) $ VU.map () $ R.toUnboxed gradQ1
 in mol{getForce = newForce, getEnergy = [energy1,energy2]} -}  
 

calcQuadraticEnergy ::  DeltaQ -> Energy -> GradMtx -> HessMtx ->  (Double,Grad)
calcQuadraticEnergy deltaQ e0 grad0 hess = 
  let energy =  e0+ termGrad + 0.5*termHess
      calcgrad = R.zipWith (+) grad0 gradInt
      grad = computeUnboxedS $ R.reshape (Z:. dim) calcgrad 
  in (energy,grad)

  where gradInt = mmultS hess deltaQ
        termGrad = sumAllS $ mmultS deltaT grad0
        termHess = sumAllS $ mmultS deltaT $ gradInt
        (Z:.dim:.1) = extent deltaQ
        deltaT   = computeUnboxedS $ reshape (Z:.1:.dim) deltaQ     