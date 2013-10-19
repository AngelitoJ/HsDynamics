{-# LANGUAGE FlexibleContexts,BangPatterns #-}

module LinearAlgebra ( 
       EigenValues
      ,EigenVectors
      ,EigenData(..)
      ,Vec(..)
      ,dimTriang
      ,dotMatrix
      ,flattenZero
      ,fourIndex2Flat
      ,identity
      ,indexDIM2toFlat
      ,indexFlat2DIM2
      ,list2ArrDIM1
      ,list2ArrDIM2
      ,mmultP
      ,mmultS
      ,mmultDIM2FlattenP
      ,mmultFlattenDIM2P
      ,mmultFlattenP
      ,mmultFlattenS
      ,scalarxMtx
      ,rowFlattenArray
      ,takeZero
      ,toTriang
      ,trace
      ,transpose2P
      ,triang2DIM2
      ,vec2Diagonal
      ,vector2Matrix
      ,zero
      ) where

-- ########## MODULE FOR IMPLENTING LinearAlgebra Tools #################


import Data.List (lookup, tails, transpose, zip5)
import Control.Applicative
import qualified Data.Vector.Unboxed as VU
import qualified Data.Map as M
import Data.Array.Repa          as R
import Data.Array.Repa.Unsafe  as R
import Control.Monad(ap,mplus)


-- ========================> TYPES and CLASSES <===================
type DIM = Int
type Indx = (Int,Int)
type Tolerance = Double
type Step = Int
type EigenValues = VU.Vector Double
type EigenVectors = Array U DIM2 Double
data EigenData = EigenData {
             eigenvals :: !EigenValues
           , eigenvec :: !EigenVectors } deriving (Show)


-- ===============> DATA TYPES <=================

newtype Vec a = Vec {runVec :: [a]} deriving Show

-- ===============> INSTANCES <==================

instance Functor Vec where
  fmap f (Vec v) = Vec $ f `fmap` v
{-
class (Functor f) => Applicative f where
pure :: a -> f a
(<*>) :: f (a -> b) -> f a -> f b
-}
instance Applicative Vec where
  pure x = Vec $ (repeat x)
  Vec fs <*> Vec xs = Vec $ Prelude.zipWith ($) fs xs
  
-- ================================== > Functions take from Repa-Examples  <==================
-- | Matrix matrix multiply.
mmultP  :: Monad m
        => Array U DIM2 Double
        -> Array U DIM2 Double
        -> m (Array U DIM2 Double)

mmultP arr brr
 = do   trr      <- transpose2P brr
        let (Z :. h1  :. _)  = extent arr
        let (Z :. _   :. w2) = extent brr
        computeP
         $ fromFunction (Z :. h1 :. w2)
         $ \ix   -> R.sumAllS
                  $ R.zipWith (*)
                        (unsafeSlice arr (Any :. (row ix) :. All))
                        (unsafeSlice trr (Any :. (col ix) :. All))
{-# NOINLINE mmultP #-}

mmultS  :: Array U DIM2 Double
        -> Array U DIM2 Double
        -> Array U DIM2 Double
mmultS arr brr
 = let trr = R.transpose brr
       (Z :. h1  :. _)  = extent arr
       (Z :. _   :. w2) = extent brr
    in computeUnboxedS
         $ fromFunction (Z :. h1 :. w2)
         $ \ix   -> R.sumAllS
                  $ R.zipWith (*)
                        (unsafeSlice arr (Any :. (row ix) :. All))
                        (unsafeSlice trr (Any :. (col ix) :. All))
{-# NOINLINE mmultS #-}

-- | Transpose a 2D matrix.
transpose2P
        :: Monad m
        => Array U DIM2 Double
        -> m (Array U DIM2 Double)

transpose2P arr
 = computeUnboxedP
 $ unsafeBackpermute new_extent swap arr
 where  swap (Z :. i :. j)      = Z :. j :. i
        new_extent              = swap (extent arr)
{-# INLINE transpose2P #-}


-- | Take the row number of a rank-2 index.
row :: DIM2 -> Int
row (Z :. r :. _) = r
{-# INLINE row #-}


-- | Take the column number of a rank-2 index.
col :: DIM2 -> Int
col (Z :. _ :. c) = c
{-# INLINE col #-}

-- Do the dot product between two matrices
dotMatrix :: Array U DIM2 Double -> Array U DIM2 Double -> Double
-- dotMatrix mtx1 mtx2 = sumAllS $ mmultS v1 v2
dotMatrix mtx1 mtx2 = sumAllS. computeUnboxedS $ fromFunction (Z:.dim^2) $
                      \(Z:. i) ->  let x = v1 ! (Z:.i)
                                       y = v2 ! (Z:.i)
                                   in  x*y

  where (Z:. dim :._) = extent mtx1
        v1 = computeUnboxedS $ reshape (Z:. dim^2) mtx1
        v2 = computeUnboxedS $ reshape (Z:. dim^2) mtx2
{-# INLINE dotMatrix #-}

        
-- ==============> SIMMETRIC MATRIX MULTIPLICATION <===============

mmultFlattenP :: Monad m => Array U DIM1 Double -> Array U DIM1 Double -> m (Array U DIM2 Double)
mmultFlattenP !arr !brr =
   do   computeUnboxedP
         $ fromFunction (Z :. dim :. dim)
         $ \(Z:. i:. j) ->  R.sumAllS
                            $ R.zipWith (*)
                              (rowFlattenArray arr i)
                              (rowFlattenArray brr j)

  where (Z:. h1) = extent arr
        dim =  dimTriang h1
{-# NOINLINE mmultFlattenP #-}

mmultFlattenS :: Array U DIM1 Double -> Array U DIM1 Double -> Array U DIM2 Double
mmultFlattenS !arr !brr =
   do   computeUnboxedS
         $ fromFunction (Z :. dim :. dim)
         $ \(Z:. i:.j ) ->   R.sumAllS
                             $ R.zipWith (*)
                               (rowFlattenArray arr i)
                               (rowFlattenArray brr j)

  where (Z:. h1) = extent arr
        dim =  dimTriang h1
{-# NOINLINE mmultFlattenS #-}   

mmultDIM2FlattenP :: Monad m => Array U DIM2 Double -> Array U DIM1 Double -> m (Array U DIM2 Double)
mmultDIM2FlattenP !arr !flat =
   do   computeUnboxedP
         $ fromFunction (Z :. dim :. dim)
         $ \ix@(Z:. i:.j) -> R.sumAllS
                             $ R.zipWith (*)
                             (unsafeSlice arr (Any :. (row ix) :. All))
                             (rowFlattenArray flat j)

  where (Z:. dim :. _) = extent arr
{-# NOINLINE mmultDIM2FlattenP #-}

mmultFlattenDIM2P :: Monad m => Array U DIM1 Double -> Array U DIM2 Double -> m (Array U DIM2 Double)
mmultFlattenDIM2P !flat !brr = do
   trr <- transpose2P brr
   computeP
     $ fromFunction (Z :. dim :. dim)
     $ \ix@(Z:. i:.j) -> R.sumAllS
                         $ R.zipWith (*)
                         (rowFlattenArray flat i)
                         (unsafeSlice trr (Any :. (col ix) :. All))

  where (Z:. dim:. _) = extent brr
{-# NOINLINE mmultFlattenDIM2P #-}

-- function to take the raw of a symmetric matrix  
rowFlattenArray :: Array U DIM1 Double -> Int -> Array U DIM1 Double
rowFlattenArray !arr !k =
  let (Z:. h1) = extent arr
      dim = dimTriang h1
  in computeUnboxedS $ fromFunction (Z:.dim)
    $  \(Z:.i) -> let x = indexDIM2toFlat dim k i  in arr ! (Z:.x)
{-# INLINE rowFlattenArray #-}

-- ===================> <=============================================
scalarxMtx :: Array U DIM1 Double -> Double -> Array U DIM1 Double
scalarxMtx !arr !s = R.computeUnboxedS $ R.unsafeTraverse arr id (\f sh -> s * f sh)
{-# INLINE scalarxMtx #-}
                   
-- ==============> FUNCTIONS FOR GENERATING REPA ARRAYS <=====================

identity :: VU.Unbox Double  => Int -> Array U DIM2 Double
identity dim =  list2ArrDIM2 dim [dij i j | i <- [1..dim], j <- [1..dim]]
  where dij u v = if u == v then 1.0 else 0.0
  
zero :: Int -> Array U DIM2 Double
zero dim = R.computeUnboxedS $ R.fromFunction (Z:. dim:. dim)
           $ (\sh -> 0.0)

flattenZero :: Int -> Array U DIM1 Double
flattenZero n = R.computeUnboxedS $ R.fromFunction (Z:. dim)
                $ (\sh -> 0.0)
  where dim = (n^2 +n) `div`2              

list2ArrDIM2 ::VU.Unbox a => Int -> [a] -> Array U DIM2 a
list2ArrDIM2 dim !list = R.fromListUnboxed (Z:. dim :. dim :: DIM2) list

list2ArrDIM1 ::VU.Unbox a => Int -> [a] -> Array U DIM1 a
list2ArrDIM1 dim !list = R.fromListUnboxed (Z:. dim :: DIM1) list
        
vec2Diagonal :: Array U DIM1 Double -> Array U DIM2 Double
vec2Diagonal vec = computeUnboxedS $ fromFunction (Z:.dim :. dim)
   (\(Z:. i:. j) -> case i == j of
                         True -> vec ! (Z:. i)
                         False -> 0.0)

  where (Z:.dim) = R.extent vec


-- ===================> Triangular Matrices <====================

-- | change the representation from Array U DIM1 a to Array U DIM2 a

takeZero :: Int -> [Int] -> [Int]
takeZero !n !xs | null xs   = [] 
                | n == 0    = [0]
                | otherwise = take n xs 
                 
triang2DIM2 :: (VU.Unbox a, Monad m) => Array U DIM1 a -> m (Array U DIM2 a)
triang2DIM2 ar = R.computeUnboxedP $ R.fromFunction (Z:.dim :. dim)
                 $ \sh@(Z:.i :.j) ->
                      let x = indexDIM2toFlat dim i j 
                      in ar ! (Z:.x)
                    
  where (Z:. len) = R.extent ar
        dim = dimTriang $ len
                
toTriang :: Monad m => Array U DIM2 Double -> m (Array U DIM1 Double)
toTriang arr = R.computeUnboxedP $ R.fromFunction (Z:. triang)
              (\(Z:. i) ->
              let (x,y) = indexFlat2DIM2 dim i
              in arr ! (Z:. x:. y) )
                 
  where (Z:.dim:. _) = extent arr                  
        triang = sum [dim,pred dim .. 1]

trace :: Array U DIM2 Double -> Double
trace arr = sumAllS  $ (R.computeUnboxedS $ R.fromFunction (Z:. dim) $
                               \(Z:.i) -> arr ! (Z:. i :. i))

  where (Z:. dim :. _) = extent arr
                               
vector2Matrix :: Int -> [a] -> [[a]]
vector2Matrix !n !list = [fmap (\i -> list !! (xs+i) ) [0..pred n] | xs <- xss]
  where xss = [n*k | k <-[0..pred n] ]
-- ============> Functions to manipulate indexes <================

indexFlat2DIM2 :: Int -> Int -> (Int,Int)
indexFlat2DIM2  !dim !i = (x,y)
  where x = f i 0
        y = g i x
        f !x !n = if x < (dim-n) then n else f (x-(dim-n)) (succ n)
        g !i !x =  i + x - (calcElem x)
        calcElem x = sum . takeZero x $ iterate pred dim


indexDIM2toFlat :: Int -> Int -> Int -> Int 
indexDIM2toFlat !dim !x !y | x <= y = calcElem x + (y-x)
                           | otherwise = calcElem y + (x-y)
  where calcElem w = sum . takeZero w $ iterate pred dim
        

sumIndexTriang:: Int -> Int
sumIndexTriang !n = (n^2 + n) `div` 2
                   
sumFirstIndex :: Int -> Int                   
sumFirstIndex !n  = n + (sumIndexTriang $ pred n) + (sum [sumIndexTriang n - k | k <- [1..pred n]])                  

totalFirstIndex :: Int -> Int -> Int
totalFirstIndex !n !firstIndex = sum [sumFirstIndex (n-m) | m <- [0..pred firstIndex]]

-- | Function to  retrieve the index of the elem n, in the flatten array (DIM1)
-- | representation containing the electronic repulsion integral (ERIs) of the form
-- | <ab|cd>
fourIndex2Flat :: Int -> [Int] -> Int
fourIndex2Flat !n !xs@[i,j,k,l] = a0 + (selectPosition n xs)
  where a0 = totalFirstIndex n i
  
selectPosition :: Int -> [Int] -> Int
selectPosition !n !xs@[i,j,k,l] =  
  case test of
     Just 1 -> l - i
     Just 2 -> (n - i) + nfun (ipred n) (ipred k) (ipred l)
     Nothing   -> (n - i) + (sumIndexTriang $ (pred n) - i) + 
                  (sum [sumIndexTriang (n-i) - m | m <- [1..(pred j - i)]])  +
                  (nfun n k l - nfun n i j)
     
  where test = first3equal `mplus` first2equal
        first3equal = if all (i==) [j,k] then Just 1 else Nothing
        first2equal = if  i == j then Just 2 else Nothing
        ipred x = iterate pred x !! (succ i)
        nfun = indexDIM2toFlat  
 


-- ================> Miscellaneous <==============================
                          
dimTriang :: Int -> Int
dimTriang a = (-1+p) `div` 2
  where p =  floor . sqrt . fromIntegral $ 1 + 8*a


     