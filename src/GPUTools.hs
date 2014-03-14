
module GPUtools where


import Data.Array.Accelerate as A
import Data.Array.Accelerate.Interpreter(run)


-- ================> <==================

main = do
  print prueba

prueba :: Array DIM1 Double
prueba = let arr1 = enumFromN (index1 100) 0
         in run $ A.map (*100) arr1
         
         
         