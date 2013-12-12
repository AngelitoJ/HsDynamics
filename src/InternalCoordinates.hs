{-# Language BangPatterns #-}


module InternalCoordinates (
                           angle
                          ,bond
                          ,calcInternals
                          ,chunks
                          ,dihedral
                          ,transform2Cart
                           ) where

import Control.Applicative             
import Control.Monad
import Control.Lens
import Control.Parallel.Strategies (parMap,rseq)
import Data.Array.Repa as R hiding ((++))
import Data.List as DL
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU

-- =======> Internal Modules <========
import CommonTypes
import LinearAlgebra


-- =============> TYPES  <=================
type VD = VU.Vector Double


-- ================> API <==========
calcInternals :: Connections -> Molecule -> Internals
calcInternals conex mol =  VU.generate dim fun
  where vu    = mol ^. getCoord . to R.toUnboxed
        vss   = chunks 3 vu
	dim   = V.length conex
	ind x = vss V.! x 
        fun x = case (conex V.! x) of
                     Bond     a b     -> bond     (ind a) (ind b) 
                     Angle    a b c   -> angle    (ind a) (ind b) (ind c)
                     Dihedral a b c d -> dihedral (ind a) (ind b) (ind c) (ind d)
                     
transform2Cart :: Connections -> Internals -> Array U DIM1 Double -> IO Grad
transform2Cart conex gradInt cart = do 
    let xss       = chunks 3 $ R.toUnboxed cart
        dim       = VU.length gradInt
        numat     = V.length xss
        wilsonMtx = calcWilson conex numat $ xss 
    wilsonT   <- transpose2P wilsonMtx 
    gradCart <- mmultP wilsonT $ R.fromUnboxed (Z:. dim :. 1) gradInt
    computeUnboxedP $ R.reshape (Z :. 3*numat) gradCart
        

calcWilson :: Connections -> Int ->  V.Vector VD -> Array U DIM2 Double
calcWilson conex numat vss = R.fromUnboxed (Z :. dimInt :.dimCart) $ VU.concat $ parMap rseq (calcDerv vss numat) xs
  where xs      = V.toList conex
        dimInt  = V.length conex
        dimCart = 3*numat
        
        
calcDerv :: V.Vector VD -> Int -> InternalCoord -> VD
calcDerv vss numat q = VU.concat $ fmap filler [0..pred numat]

  where ind x    = vss V.! x 
        zeros    = VU.generate 3 $ const 0
        look     = fromMaybe (error "There was a problem calculating the Wilson matrix") . flip DL.lookup dervs
        filler x = if x`elem` ixs then look x else zeros 
        (ixs,dervs) = case q of
                        Bond     a b     -> let xs = [a,b]
                                            in (xs, zip xs $ derv_bond (ind a) (ind b)) 
                        Angle    a b c   -> let xs = [a,b,c]
                                            in (xs,zip xs $ derv_angle (ind a) (ind b) (ind c))                                           
                        Dihedral a b c d -> let xs = [a,b,c,d]
                                            in (xs,zip xs $ derv_dih (ind a) (ind b) (ind c) (ind d)) 
                     
                       
-- =============> <==================
vecdot :: VD -> VD -> Double
vecdot v1 v2 = VU.sum $ VU.zipWith (*) v1 v2

vecnorm :: VD -> Double
vecnorm v = sqrt $ vecdot v v

vecCross :: VD -> VD -> VD  
vecCross !v1 !v2 = VU.fromList $ fmap permute [(1,2),(2,0),(0,1)]
  where permute (i,j) = (v1 VU.! i)*(v2 VU.! j) - (v2 VU.! i)*(v1 VU.! j)

vecsum :: VD -> VD -> VD
vecsum !v1 !v2 = VU.zipWith (+) v1 v2        
  
vecsub :: VD -> VD -> VD
vecsub !v1 !v2 = VU.zipWith (-) v1 v2        
        
dif2 :: VD -> VD -> VD -> Double
dif2 !p1 !p2 !p3 =  (p1 `vecsub` p2) `vecdot` (p3 `vecsub` p2)


vecScalar :: VD -> Double -> VD
vecScalar vd s = VU.map (*s) vd

-- =============> Bonds, Angles and Dihedrals <============                 
bond :: VD -> VD -> Double
bond !p1 !p2 = vecnorm $ p1 `vecsub` p2

angle :: VD -> VD -> VD -> Double
angle !p1 !p2 !p3 = acos $ arg
  where arg = (dif2 p1 p2 p3) / ((p1 `bond` p2) * (p2 `bond` p3))

dihedral :: VD -> VD -> VD -> VD -> Double
dihedral !p1 !p2 !p3 !p4  =
  let [xba,xca,xcb,xdb] = DL.zipWith (vecsub) [p2,p3,p3,p4] [p1,p1,p2,p2]
      [w1,w2] = DL.zipWith vecCross [xba,xcb] [xca,xdb]
      [n1,n2] = fmap vecnorm [w1,w2]
      teta = acos $ ((w1 `vecdot` w2) / (n1*n2))
  in if  (signum  $ w2 `vecdot` xba) > 0.0  then teta else -teta

-- ===========>  Internal Coordinates derivatives <======================

derv_bond :: VD -> VD -> [VD]
derv_bond a b = [da,mda] 
  where r = bond a b
        da  = VU.map (/r) $ vecsub a b
        mda = VU.map negate da
        

derv_angle :: VD -> VD -> VD -> [VD]      
derv_angle a b c = [vecsum a1 a2,vecsub b1 b2, vecsum c1 c2]
  where [r1,r2] = DL.zipWith bond [a,b] [b,c]
        teta    = angle a b c
        [d1,d2] = DL.zipWith vecsub [b,a] [c,b]
        d3      = let b2a = vecsub (VU.map (*2) b) a in vecsub b2a c
        f1      = recip $ r1*r2*(sin teta)               
        f2      = recip $ tan teta
        a1      = VU.map (*f1) d1
        a2      = VU.map (*(f2/r1^2)) d2
        b1      = let x = VU.map (negate . (*(f2/r1^2))) d2 
                      y = VU.map (/r2^2) d1
                  in vecsum x y
        b2      = VU.map (*f1) d3 
        c1      = VU.map (negate . (*f1)) d2
        c2      = VU.map (negate . (*(f2/r2^2))) d1        

derv_dih :: VD -> VD -> VD -> VD -> [VD]
derv_dih a b c d = [da,db,dc,dd]
  where r2         = bond b c
        [v1,v2,v3] = DL.zipWith vecsub [a,c,c] [b,b,d]
        [uxw,vxw]  = DL.zipWith vecCross [v1,v2] [v2,v3]
        om1        = vecdot uxw uxw
        om2        = vecdot vxw vxw
        v1_v2      = vecdot v1 v2
        v3_v2      = vecdot v3 v2
        p1         = v1_v2/r2^2 - 1
        p2         = v3_v2/r2^2
        p3         = p2 -1
        p4         = p1 +1
        
        da = VU.map (*(r2/om1)) uxw
        db = vecsub (vecScalar da p1) (vecScalar dd p2)
        dc = vecsub (vecScalar dd p3) (vecScalar da p4)
        dd = VU.map (negate . (*(r2/om2))) vxw
      
-- ================> Utilities <=====================    
  
chunks :: Int -> VD ->  V.Vector VD
chunks n v = v1 V.++  v2 
  where v1 = V.map (\x -> VU.backpermute v $ VU.generate n $ \y -> n*x +y ) ixs
        v2 = if r == 0 then V.empty else V.singleton $ VU.backpermute v $ VU.generate r $ \i -> (n*q) + i
        len = VU.length v
        (q,r) = quotRem len n
        ixs = V.generate q id
        

        
            
  


  