

module Prueba where

import Gaussian 

-- ===========================> <==========


main =  do
  eitherResult <- parseGaussianCheckpoint "watersto3G.fchk"
  case eitherResult of
       Left x -> print x
       Right w -> mapM_ fun w 
       
       
  where fun = \r -> case r of
                       IGauBlock label _n xs -> print label
                       RGauBlock label _n xs -> print label
                       otherwise -> return ()
                       
                       
                       
        
  
  