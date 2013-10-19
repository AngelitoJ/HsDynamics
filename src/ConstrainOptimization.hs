
module ConstrainOptimization where

import Control.Applicative
import Data.Array.Repa as R
import qualified Data.List as DL
import Numeric.LinearAlgebra as NLA((|>),(<>),(<.>),(@>), Matrix,Vector,add,buildMatrix,
                                   eigenvalues,eigSH',fromColumns,fromList,fromLists,fromRows,inv,
                                   rows,scale,sub,toList,trans)

-- Internal Modules                                   
import CommonTypes
import Dynamics
import APIparser
import Constants (auN)
                                   
-- ==============> <================
type ExternalGrad = Array U DIM1 Double
type Hessian = Matrix Double
type Gradient = Vector Double
type Position = Vector Double
type Step = Int
                                   
-- ===============>   <============

driverConstrainedOptimizer :: Job -> String -> Molecule ->  (Int,Int) -> Double -> IO ()
driverConstrainedOptimizer job project mol (n,m) norm = do
--   let gext = computeUnboxedS . R.map (*(-1)) $ calcFext (n,m) (norm/auN) mol
  let gext = calcFext (n,m) (norm/auN) mol
      hess = buildMatrix (dim) (dim) (\(i,j) -> if i==j then 1 else 0)
      dij x y = if x == y then 1.0 else 0.0
      (Z:. dim) =  extent . getForce $ mol
  result <- mainLoop job project mol hess gext 0
  putStrLn "successful convergence!!"
                          
mainLoop :: Job -> String -> Molecule -> Hessian -> ExternalGrad -> Step -> IO Molecule
mainLoop job project mol0 hess gext step = do
   let (mol1,dx)  = updateGeometry mol0 hess gext 
   newmol <- interactWith job project mol1
   let [oldGrad,newGrad] = fmap  getForce [mol1,newmol]
       deltaGrad = NLA.fromList . R.toList . computeUnboxedS $ R.zipWith (-) newGrad oldGrad
       newHess = updateHess dx deltaGrad hess
       gradT = totalGrad newmol gext
       norm = sqrt . R.sumAllS . computeUnboxedS $ R.zipWith (*) gradT gradT
   putStrLn $ "Step number: " DL.++ (show step)
   putStrLn $ "Norm total Gradient au: " DL.++ (show norm)
   printMol newmol $ "total Gradient (au): " DL.++ (show norm)
   if converge norm  then return newmol
                     else mainLoop job project newmol newHess gext (succ step)
  
totalGrad :: Molecule -> ExternalGrad -> Array U DIM1 Double
totalGrad mol gext = R.computeUnboxedS . R.zipWith (+) gext $ getForce mol  

updateGeometry :: Molecule -> Hessian -> ExternalGrad -> (Molecule,Position)                                             
updateGeometry mol hess gext = (mol{getCoord = newGeom},dx)
  where newGeom = R.computeUnboxedS $ R.zipWith (+) oldgeom dxRepa
        oldgeom = getCoord mol
        dxRepa =  R.fromListUnboxed (extent oldgeom) . NLA.toList $ dx
        dx =  scaleDX $ newtonRaphson grad hess
        grad = R.toList $ totalGrad mol gext
                      
newtonRaphson :: [Double] -> Hessian -> Position
newtonRaphson gs hs =  negateInv <> (fromList gs)
  where  negateInv = scale (-1.0) hInverse
         hInverse = inverseNonNegative hs

inverseNonNegative :: Hessian -> Hessian
inverseNonNegative hess = inv $ hess `add` (scale l unit)
  where  (eigVals,_) = eigSH' hess
         vals = fmap negate . NLA.toList  $ eigVals
         dim = rows hess
         unit = buildMatrix dim dim (\(i,j)-> if i ==j then 1 else 0)
         l = maximum  (0 : vals)
         

calcFext :: (Int,Int) -> Double -> Molecule -> Array U DIM1 Double
calcFext (n,m) normFext mol = R.fromListUnboxed sh xs
  where  xs = concatMap (sparseList m n deltaF) [0..numat]
         [v1,v2] = fmap (\j -> let i = 3*(pred j) in fmap (coords !) [(Z:.i),(Z:.i+1),(Z:.i+2)]) [m,n]
         deltaF = fmap (*normFext) . normalize $ DL.zipWith (-) v2 v1
         sh@(Z:. dim) = extent coords
         coords = getCoord mol
         numat = pred (dim `div` 3)
       
         
sparseList :: Int -> Int -> [Double] -> Int -> [Double]
sparseList m n xs i | i == (pred m) = xs
                    | i == (pred n) = fmap negate xs
                    | otherwise = take 3 . repeat $ 0.0

        
updateHess :: Position -> Gradient -> Hessian -> Hessian
updateHess dx dgrad hess =  hess `add` termGrad `sub` termHess
  where termGrad = scale (recip $ dgrad <.> dx) $ (toMatrixColumn  dgrad) <> (toMatrixRow dgrad)
        termHess = scale (recip dot) $ mtx1 <> mtx2
        dot = dx <.> (hess <> dx)
        mtx1 = hess <> (toMatrixColumn dx)
        mtx2 = trans mtx1

         
converge :: Double -> Bool
converge norm = if norm < 1.0e-4 then True else False 

scaleDX :: Position -> Position
scaleDX dx = if norm < 0.5 then dx 
                           else scale (0.5/norm) dx
  where norm = sqrt $ dx <.> dx

-- =======================================================
toMatrixRow :: Vector Double -> Matrix Double
toMatrixRow = fromRows . pure 

toMatrixColumn :: Vector Double -> Matrix Double
toMatrixColumn = fromColumns . pure
 
normalize :: [Double] -> [Double]
normalize xs = fmap (*(recip norm)) xs
  where norm = sqrt . sum $ DL.zipWith (*) xs xs        
                                  

diagonal2Mtx :: Int -> [Double] -> Matrix Double
diagonal2Mtx numat xs = buildMatrix (3*numat) (3*numat) fun 
  where dim = pred $ 3*numat
        fun (i,j) = if i<=j then fun2 i j else fun2 j i
        fun2 x y = xs !! (calcElem x + (y-x))
        calcElem x = sum . takeZero x $ iterate pred (3*numat)       
        
takeZero :: Int -> [Int] -> [Int]
takeZero n xs | null xs   = [] 
              | n == 0    = [0]
              | otherwise = take n xs       