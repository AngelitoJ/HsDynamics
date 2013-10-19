
-- HsDynamics: 

-- @2013 Felipe Zapata, Alessio Valentini, Angel Alvarez from The ResMol Group  

module Main where


import Data.List (zipWith3)
import Data.Complex
import Data.Maybe ( fromMaybe )
import Control.Concurrent
import Control.Concurrent.Async
import Control.Monad ((<=<),liftM)
import Control.Monad.Trans.Either
import System.Environment ( getArgs )
import System.FilePath
import System.Cmd ( system )
import System.Console.GetOpt
import Text.ParserCombinators.Parsec
import Text.Printf


-- Cabal imports
import Data.Version (showVersion)
import Distribution.Version
import Paths_HsDynamics as HsDynamics

-- internal imports
import APIparser
import BasisParser
import CommonTypes
import Constants  
import ConstrainOptimization
import Dynamics
import InitialConditions
import Gaussian
import GenericParser
import Molcas
import OptsCheck
import Tasks
import Tully



program = "Molecular Dynamics"
authors = "@2013  Felipe Zapata, Alessio Valentini, Angel Alvarez"

-- default options
defaultOptions    = Options
 { optDump        = False
 , optModules     = [("constrained",processConstrained),("dynamics",processDynamics),("externalForces",processExternalForces),("molcas",processMolcas),("generic",processGenericFiles),("basis",processBasisFiles)]
 , optMode        = Nothing
 , optVerbose     = False
 , optShowVersion = False
 , optOutput      = Nothing
 , optDataDir     = Nothing
 , optInput       = []
 , optTemperature = Nothing
 }

-- currently supported options
acceptedOptions :: [OptsPolicy]
acceptedOptions =
 [ 
   Option ['h','?'] ["help"]    (NoArg  ( check_help           ))                "Show this help message."
 , Option ['v']     ["verbose"] (NoArg  ( check_verbosity      ))                "Verbose run on stderr"
 , Option ['V']     ["Version"] (NoArg  ( check_version        ))                "Show version number"
 , Option ['D']     ["0"]       (ReqArg ( check_data_dir       ) "Dir")          "Directory where files are located"
 , Option ['m']     ["mode"]    (ReqArg ( check_operation_mode ) "Mode")         "Mode of Operation"
 , Option []        ["dump"]    (NoArg  ( check_dump_options   ))                "Force args cmdline dump"
 , Option ['t']     ["temperature"] (ReqArg (check_temperature) "Temperature")   "Temperature of the simulation"
  ]
--    Option ['e']     ["error"]   (NoArg (\ _opts -> return $ Left "forced error on args detected!"))  "Force args checking error"
--  , Option ['i']     ["input"]   (OptArg (\f opts -> check_input_file f opts) "FILE")             "Input file"


main :: IO ()
main = do
    args   <- getArgs
    cores  <- getNumCapabilities
    progHeader cores
    result <- runEitherT $ progOpts args defaultOptions acceptedOptions
    either somethingIsWrong doSomeStuff result


somethingIsWrong :: String -> IO ()    
somethingIsWrong msg = do
             putStrLn $ "\nError: " ++ msg ++ "\n"
             putStrLn $ usageInfo header acceptedOptions

doSomeStuff :: Options -> IO ()
doSomeStuff optsR@Options { optMode = mode } = do
    case mode of
         Nothing -> printFiles optsR
         Just fun -> fun optsR

-- Keep calm and curry on, we are the good guys....
progHeader :: Int -> IO ()
progHeader c = 
    putStrLn $ program ++ " V:" ++ currVersion ++ " " ++ authors ++ "\n\t" ++ show(c) ++ " processor " ++ (core2string c) ++ " detected."
    where
        currVersion :: String
        currVersion = showVersion HsDynamics.version
        core2string :: Int -> String
        core2string c = case c > 1 of
                             True -> "cores"
                             False -> "core"

header :: String
header = "Usage: Options [OPTION...] files..."

-- | "Efects for dummies", this functions has no purpouses other than printng args
printFiles :: Options -> IO ()
printFiles opts@Options { optInput = files, optDataDir = datadir } = do
    putStrLn $ "Processing args with options:\n" ++ (show opts) ++ "\n\n"
    mapM_ printargs filepaths 
    where
            dir = fromMaybe "" datadir
            filepaths = zipWith (combine) (cycle [dir]) files
            printargs :: String -> IO ()
            printargs path = putStrLn $ "Processing path: " ++ path ++ "..."

processGenericFiles :: Options -> IO ()
processGenericFiles opts@Options { optInput = files, optDataDir = datadir } = do
    mapM_ processGenericFile filepaths 
    where
            dir = fromMaybe "" datadir
            filepaths = zipWith (combine) (cycle [dir]) files

processBasisFiles :: Options -> IO ()
processBasisFiles opts@Options { optInput = files, optDataDir = datadir } = do
    mapM_ processBasisFile filepaths 
    where
            dir = fromMaybe "" datadir
            filepaths = zipWith (combine) (cycle [dir]) files


-- =============> Drivers to run the molecular dynamics simulations in Molcas <==============
processMolcas :: Options -> IO ()
processMolcas opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[xyz,molcas,input] =  optInput opts
  InitialDynamics st time dt externalForce anchor project <- parseInputDynamics input
  initialMol <- initializeMolcasOntheFly xyz st temp
  let numat = length $ getAtoms initialMol
      [auTime,audt] = fmap (/au_time) [time,dt]
      thermo = initializeThermo numat temp
      step = 10
      job = Molcas $ project 
      aMatrix = initialAMTX initialMol
  mol <- interactWith job project initialMol
  constantForceDynamics mol job thermo temp auTime audt anchor externalForce aMatrix step
  
  
  
            
-- =============> Drivers to run the molecular dynamics simulations in Gaussian <==============

-- | on the fly molecular dynamics with applied external forces
processExternalForces :: Options -> IO ()
processExternalForces opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      [input,fchk,out] = optInput opts 
  InitialDynamics st time dt externalForce anchor theoryLevels <- parseInputDynamics input
  mol <- (updateMultiStates out) <=< (initializeSystemOnTheFly fchk st) $ temp
  let numat = length $ getAtoms mol
      thermo = initializeThermo numat temp
      [auTime,audt] = fmap (/au_time) [time,dt]
      aMatrix = initialAMTX mol 
      step = 1
      job = Gaussian theoryLevels
  constantForceDynamics mol job thermo temp auTime audt anchor externalForce aMatrix step
  
constantForceDynamics ::  Molecule -> Job -> Thermo -> Temperature -> Time -> DT-> Anchor -> Double -> MatrixCmplx -> Int -> IO ()
constantForceDynamics mol job thermo temp time dt anchor externalForce aMatrix step = do
   if time < 0.0 then return ()
                 else do 
                  let es = concatMap (printf "%.6f  ") . head . getEnergy $ mol
                  printMol mol es
                  printData mol step
                  (newMol,newThermo) <- dynamicExternalForces mol dt temp thermo job "TullyExternalForces" anchor externalForce
                  (tullyMol,newAmatrix) <- tullyDriver dt aMatrix step newMol
                  printGnuplot newAmatrix tullyMol
                  constantForceDynamics tullyMol job newThermo temp (time-dt) dt anchor externalForce newAmatrix (succ step)
  
  
tullyDriver ::  DT -> MatrixCmplx -> Int ->  Molecule -> IO (Molecule,MatrixCmplx)
tullyDriver dt aMatrix step mol =
  if (length $ getCoeffCI mol) /= 3  -- At least 3 set of CI coefficients are required in oder to initialize the Tully
     then return (mol,aMatrix) 
     else tullyHS dt aMatrix step mol
                
 
    
  
-- Felipe hizo que la mierda se manifestara en esta función tan horrible!!!
processDynamics :: Options -> IO ()
processDynamics = undefined
-- processDynamics  opts@Options { optInput = (fileHess2:fileGrad2:fileHess1:fileGrad1:filestate2 :filestate1:input:_)} =  
--   do [n,temp,time,dt] <- ((fmap readDouble) . words ) `fmap` (readFile input)
--      let numat  = floor n
--          dimInt = 3*n-6
--          [auTime,auDt] = fmap (/au_time) [time,dt]
--          thermo = initializeThermo numat temp  
--          initialRC = initializeRC numat
--      [st1,st2] <- mapM (takeInformation <=< parseGaussianCheckpoint) [filestate1,filestate2]
--      let ParseInformation maybeXs maybeEs1 _maybeGrad1 _maybeHess1  maybeMasses _maybeLabels1 = st1
--          ParseInformation _maybeXs maybeEs2 _maybeGrad2 _maybeHess2 _maybeMasses2 _maybeLabels2 = st2
--          [coord,[energy1],[energy2],masses] = fmap (fromMaybe (error "something wrong reading fchk files")) [maybeXs,maybeEs1,maybeEs2,maybeMasses]
--          aumasses = fmap (*amu) masses
--      conex <- parseConnections "internas.dat"
--      gs <- mapM readArrayDIM1FIle [fileGrad1,fileGrad2]
--      hs <- mapM readArrayDIM2FIle [fileHess1,fileHess2]
--      let derivatives = zipWith3 EnergyDerivatives [energy1,energy2] gs hs
--      mol <- initialConditions coord aumasses derivatives conex temp     
--      print initialRC 
--      print "initial conditions ready"
-- --      rc <- driverNoseHoover mol auTime auDt temp thermo Quadratic "quadratic" initialRC
--   where readDouble x = read x :: Double


driverNoseHoover :: Molecule -> Time -> DT -> Temperature -> Thermo -> Job -> String -> ReactionCoordinate -> IO ReactionCoordinate
driverNoseHoover mol time dt temp thermo job project rc =
  if time < 0.0 then return rc
                else do 
                    (newMol,newThermo) <- dynamicNoseHoover mol dt temp thermo job project
                    let newRC = updateReactionCoordinate rc newMol
                    driverNoseHoover newMol (time-dt) dt temp newThermo job project newRC
 

updateReactionCoordinate :: ReactionCoordinate -> Molecule -> ReactionCoordinate
updateReactionCoordinate = undefined



processConstrained :: Options -> IO ()
processConstrained = undefined
-- processConstrained opts@Options { optInput = (file:_)} =
--   do  let n1 = pred 16
--           n2 = pred 23
--           norm = 5.0e-9 -- force in nanoNewton 
--       mol <- parseJob defaultMol{getAtoms=labels} job file
--       r <- parseGaussianCheckpoint file          
--       driverConstrainedOptimizer job "azoCis_5nN" mol (16,23) norm   
--            
--   where job = (Gaussian theory)
--         theory = "cam-b3lyp/6-311+g(d,p)"          
--         labels = ["c","c","c","c","c","c","n","n","c","c","c","c","c","c","H","H","H","H","H","H","H","H","H","H"]
--         numat = length labels
--         fun x = case x of
--                      RGauBlock label _n _xs -> if (words label) == ["Cartesian","Force","Constants"] then True else False
--                      otherwise -> False

 

