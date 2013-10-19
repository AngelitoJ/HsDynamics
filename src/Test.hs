
-- HsDynamics, a haskell ab-initio Molecular Dynamics program
-- Angel Alvarez, Felipe Zapata 
-- @ 2013 The Resmol Group


module Main where

import System.Environment
import System.IO
import Gaussian
import Molcas
import MolcasShow

-- For a complete description of this method please refer to:
--
--

main :: IO ()
main = do
    args <- getArgs
    case args of
         [file] -> do
             result <- parseMolcasOutput file
             case result of
                  Left err  -> print err
                  Right xs  -> mapM_ (putStrLn.show) (id xs)
         _ -> do
             progName <- getProgName
             hPutStrLn stderr $ "Usage: " ++ " filename"

myfilter l= filter isSomeModule modules
    where
        isSomeModule = (isModule "alaska") `orModule` (isModule "rasscf")
        (ModuleAuto  _ _ modules) =  head $ filter isAuto l
-- myfilter = filter (isModule "alaska" `orModule` isModule "rasscf")