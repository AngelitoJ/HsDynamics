
-- HsDynamics: 

-- @2013 Felipe Zapata, Alessio Valentini, Angel Alvarez from The ResMol Group  

module Main where


import Data.List (zipWith3)
import Data.Complex
import Data.Maybe ( fromMaybe )
import Control.Concurrent
import Control.Concurrent.Async
import Control.Lens ((^.),to)
import Control.Monad ((<=<),liftM)
import Control.Monad.Trans.Either
import System.Environment ( getArgs )
import System.FilePath
import System.Cmd ( system )
import System.Console.GetOpt
import Text.Printf

-- Cabal imports
import Data.Version (showVersion)
import Distribution.Version
import Paths_HsDynamics as HsDynamics

-- internal imports
import APIparser
import CommonTypes
import Constants  
import ConstrainOptimization
import Dynamics
import InitialConditions
import InternalCoordinates
import Gaussian
import Molcas
import OptsCheck
import Tasks
import TinkerQMMM
import Tully



program = "Molecular Dynamics"
authors = "@2013  Felipe Zapata, Alessio Valentini, Angel Alvarez"

-- default options
defaultOptions    = Options
 { optDump        = False
 , optModules     = [("constrained",processConstrained),("externalForces",processExternalForces),("molcas",processMolcas),("palmeiro",processPalmeiro),("molcasTinker",processMolcasTinker),("molcasZeroVel",processMolcasZeroVelocity),("prueba",processPrueba)]
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
        core2string c = if c > 1 then  "cores"
                                 else  "core"

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

-- =============> Drivers to run the molecular dynamics simulations in Molcas <==============

processPrueba :: Options -> IO ()
processPrueba opts =   do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[xyz,ctl,input] =  optInput opts
  ctl <- getSuffixFile "." ".ctl"
  conex <- parserFileInternasCtl ctl
  initData <- parseFileInput parseInput input
  let getter   = (initData ^.)
  initialMol  <- initializeMolcasOntheFly xyz (getter getInitialState) temp
  let ints = calcInternals conex initialMol
  print ints


processMolcas :: Options -> IO ()
processMolcas opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[xyz,molcas,input] =  optInput opts
  initData <- parseFileInput parseInput input
  let getter   = (initData ^.)
      project  = getter getProject
  initialMol  <- initializeMolcasOntheFly xyz (getter getInitialState) temp
  molcasInput <- parseMolcasInputFile molcas
  let numat = initialMol ^. getAtoms . to length
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      thermo = initializeThermo numat temp
      step = 1
      job = Molcas molcasInput
      aMatrix = initialAMTX initialMol
  mol <- interactWith job project initialMol
  constantForceDynamics mol job thermo temp auTime audt (getter getForceAnchor) (getter getExtForceMod) aMatrix step project
  

processMolcasZeroVelocity :: Options -> IO ()
processMolcasZeroVelocity opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[xyz,molcas,input] =  optInput opts
  initData <- parseFileInput parseInput input
  let getter   = (initData ^.)
      project  = getter getProject
  initialMol  <- initializeMolcasZeroVel xyz (getter getInitialState) temp
  molcasInput <- parseMolcasInputFile molcas
  let numat = initialMol ^. getAtoms . to length
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      thermo = initializeThermo numat temp
      step = 1
      job = Molcas molcasInput
      aMatrix = initialAMTX initialMol
  mol <- interactWith job project initialMol
  constantForceDynamics mol job thermo temp auTime audt (getter getForceAnchor) (getter getExtForceMod) aMatrix step project  
  
  
processMolcasTinker :: Options -> IO ()
processMolcasTinker opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[tinkerKey,tinkerXYZ,molcasInput,input] =  optInput opts      
  initData <- parseFileInput parseInput input
  let getter = (initData ^.)
      project  = getter getProject
  tinkerQMMM <- parserXYZFile tinkerXYZ
  atomsQM    <- parserKeyFile tinkerKey
  initialMol <- initializeMolcasTinker molcasInput (getter getInitialState) temp $ length atomsQM
  let numat = initialMol ^. getAtoms . to length
      thermo = initializeThermo numat temp
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      step = 1
      tinkerCommand = "prueba"
      job = MolcasTinker atomsQM tinkerCommand
  print initialMol
--   driverMolcasTinker initialMol audt auTime thermo job project step
   
driverMolcasTinker :: Molecule -> DT -> Temperature -> Thermo -> Job -> String  -> Step -> IO ()  
driverMolcasTinker mol dt t thermo job project step = 
  if t <0 then return ()
          else do 
            let es = concatMap (printf "%.6f  ") $ mol ^. getEnergy . to head
            printMol mol es
            printData mol step
            (newMol,newThermo) <- dynamicNoseHoover mol dt t thermo job project
            driverMolcasTinker newMol dt (t-dt) newThermo job project (succ step)
    
-- =============> Drivers to call Palmeiro Interpolator <======================================

-- | Molecular Dynamics using interpolated PES 
processPalmeiro :: Options -> IO ()
processPalmeiro opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[xyz,input] =  optInput opts         
  initData <- parseFileInput parseInput input
  let getter = (initData ^.)
  initialMol  <- initializeMolcasOntheFly xyz (getter getInitialState) temp
  ctl         <- getSuffixFile "." ".ctl"
  conex       <- parserFileInternasCtl ctl
  let numat = initialMol ^. getAtoms . to length
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      thermo        = initializeThermo numat temp
      job           = Palmeiro conex ["/S0","/S1"]
  palmeiroLoop initialMol audt temp thermo job "" 1
  
palmeiroLoop :: Molecule -> DT -> Temperature -> Thermo -> Job -> String -> Step -> IO ()
palmeiroLoop  mol dt t thermo job project step = do
    print $ "Step: " ++ show step
    if t > 0 then do
                  (newMol,newThermo) <- dynamicNoseHoover mol dt t thermo job project  
                  printMol  newMol ""
                  printData newMol step
                  return ()
--                   palmeiroLoop newMol dt (t - dt) newThermo job project (succ step)
             else return ()

-- =============> Drivers to run the molecular dynamics simulations in Gaussian <==============

-- | on the fly molecular dynamics with applied external forces
processExternalForces :: Options -> IO ()
processExternalForces opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      [input,fchk,out] = optInput opts 
  initData <- parseFileInput parseInput input
  let getter = (initData ^.)
  mol <- (updateMultiStates out) <=< (initializeSystemOnTheFly fchk $ getter getInitialState) $ temp
  let numat         = mol ^. getAtoms . to length
      thermo        = initializeThermo numat temp
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      aMatrix       = initialAMTX mol 
      step          = 1
      theoryLevels  = getter getTheory
      basis         = getter getBasis
      job = Gaussian (theoryLevels,basis)
  constantForceDynamics mol job thermo temp auTime audt (getter getForceAnchor) (getter getExtForceMod) aMatrix step "TullyExternalForces"
  
constantForceDynamics ::  Molecule -> Job -> Thermo -> Temperature -> Time -> DT-> Anchor -> Double -> MatrixCmplx -> Int -> Project -> IO ()
constantForceDynamics mol job thermo temp time dt anchor externalForce aMatrix step project = do
   if time < 0.0 then return ()
                 else do 
                  let es = concatMap (printf "%.6f  ") $ mol ^. getEnergy . to head
                  printMol mol es
                  printData mol step
                  (newMol,newThermo) <- dynamicExternalForces mol dt temp thermo job project anchor externalForce
                  (tullyMol,newAmatrix) <- tullyDriver dt aMatrix step newMol
                  printGnuplot newAmatrix tullyMol
                  let [oldRoot,newRoot] = (^.getElecSt) `fmap` [newMol,tullyMol]
                      newJob            =  if oldRoot == newRoot then job else updateNewJobInput job tullyMol
                  constantForceDynamics tullyMol newJob newThermo temp (time-dt) dt anchor externalForce newAmatrix (succ step) project
  
  
tullyDriver ::  DT -> MatrixCmplx -> Int ->  Molecule -> IO (Molecule,MatrixCmplx)
tullyDriver dt aMatrix step mol =
  if (mol^. getCoeffCI . to length) /= 3  -- At least 3 set of CI coefficients are required in oder to initialize the Tully
     then return (mol,aMatrix) 
     else tullyHS dt aMatrix step mol




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

 

