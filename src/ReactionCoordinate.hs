{-# Language BangPatterns #-}


module ReactionCoordinate where

import Control.Arrow ((&&&))
import qualified Data.List as DL
import qualified Data.Map as M
import qualified Data.Vector as V
import Data.Vector.Unboxed as VU

-- =============> Internal Imports <============
import CommonTypes

-- ================> data types <==============
type Vector3D = VU.Vector Double


-- ================> <===============
        
calculateInternals :: V.Vector InternalCoord -> VU.Vector Double -> VU.Vector Double
calculateInternals !conex !cart = VU.unfoldr genInternal conex 
  where genInternal !ints = let q = calculateQ (V.head ints) cart
                                t = V.tail ints
                            in  if V.null ints then Nothing else Just (q,t)

calculateInternalsFC :: Internals -> V.Vector InternalCoord -> VU.Vector Double -> VU.Vector Double
calculateInternalsFC !fc !conex !cart = VU.unfoldr genInternal (fc,conex) 
  where genInternal (!fcq,!ints) = let hi  = V.head ints
                                       hfc = VU.head fcq
                                       q = calculateQFC hi hfc cart
                                       ti = V.tail ints
                                       tfc = VU.tail fcq
                                   in  if V.null ints then Nothing else Just (q,(tfc,ti))
                            
                            
calculateQ :: InternalCoord -> VU.Vector Double -> Double
calculateQ i cart =
  case i of
       Bond  i j        -> bond     (getAtom i) (getAtom j)
       Angle i j k      -> angle    (getAtom i) (getAtom j) (getAtom k)
       Dihedral i j k l -> dihedral (getAtom i) (getAtom j) (getAtom k) (getAtom l)
                         
  where getAtom x = VU.unsafeSlice (3*x) 3 cart 

calculateQFC :: InternalCoord ->  Double -> VU.Vector Double -> Double
calculateQFC !i !fcq cart =
  case i of
       Bond  i j        -> bond       (getAtom i) (getAtom j)
       Angle i j k      -> angle      (getAtom i) (getAtom j) (getAtom k)
       Dihedral i j k l -> dihedralFC (getAtom i) (getAtom j) (getAtom k) (getAtom l) fcq
                         
  where getAtom x = VU.unsafeSlice (3*x) 3 cart 
    
 
 
-- =================> BOND, ANGLE AND DIHEDRAL <=========
        

vecScalar :: Double -> Vector3D -> Vector3D
vecScalar s v = VU.map (*s) v
        
vecCross :: Vector3D -> Vector3D -> Vector3D
vecCross v1 v2 = VU.fromList $ fmap permute [(1,2),(2,0),(0,1)]
  where permute (i,j) = (v1 ! i)*(v2 ! j) - (v2 ! i)*(v1 ! j)
--         
vecdot :: Vector3D -> Vector3D -> Double
vecdot v1 v2 = VU.sum $ VU.zipWith (*) v1 v2

vecnorm :: Vector3D -> Double
vecnorm v1 = sqrt $ vecdot v1 v1

vecsub :: Vector3D -> Vector3D -> Vector3D
vecsub v1 v2 = VU.zipWith (-) v1 v2

dif2 :: Vector3D -> Vector3D -> Vector3D -> Double
dif2 p1 p2 p3 =  (p1 `vecsub` p2) `vecdot` (p3 `vecsub` p2)
        
    
bond :: Vector3D -> Vector3D -> Double
bond p1 p2 = vecnorm $ p1 `vecsub` p2

angle :: Vector3D -> Vector3D -> Vector3D -> Double
angle p1 p2 p3 = (180.0/pi*) . acos $ arg
  where arg = (dif2 p1 p2 p3) / ((p1 `bond` p2) * (p2 `bond` p3))
 
dihedral :: Vector3D -> Vector3D -> Vector3D -> Vector3D -> Double
dihedral p1 p2 p3 p4  =
  let [xba,xca,xcb,xdb] = DL.zipWith (vecsub) [p2,p3,p3,p4] [p1,p1,p2,p2]
      [w1,w2] = DL.zipWith vecCross [xba,xcb] [xca,xdb]
      [n1,n2] = fmap vecnorm [w1,w2]
      teta = (180.0/pi*) . acos $ ((w1 `vecdot` w2) / (n1*n2))     
  in if  (signum  $ w2 `vecdot` xba) > 0.0  then teta else -teta

dihedralFC :: Vector3D -> Vector3D -> Vector3D -> Vector3D -> Double ->  Double
dihedralFC p1 p2 p3 p4 fc =
  let [xba,xca,xcb,xdb] = DL.zipWith (vecsub) [p2,p3,p3,p4] [p1,p1,p2,p2]
      [w1,w2] = DL.zipWith vecCross [xba,xcb] [xca,xdb]
      [n1,n2] = fmap vecnorm [w1,w2]
      val = (180.0/pi*) . acos $ ((w1 `vecdot` w2) / (n1*n2))
      teta = if (signum  $ w2 `vecdot` xba) > 0.0  then val else -val
      angles = [teta,teta-360,teta+360]
      xs = DL.zip (fmap (\x -> abs $ x-fc) angles) [0..]
      ix = DL.head . fmap snd  $ DL.sort xs
  in  angles !! ix
      
      
                
                