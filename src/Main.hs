{-# Language RankNTypes #-}
-- HsDynamics: 

-- @2013 Felipe Zapata, Alessio Valentini, Angel Alvarez from The ResMol Group  

module Main where


import Data.List (zipWith3)
import Data.Complex
import Data.Maybe ( fromMaybe )
import Control.Concurrent
import Control.Concurrent.Async
import Control.Lens ((.~),(^.),(&),Getting(..),to)
import Control.Monad ((<=<),liftM,zipWithM_)
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
import Logger 
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
 , optModules     = [("gaussTully",processGauss), ("molcasTully",processMolcas),                                           
                     ("verletGaussian",processVerletGaussian),("verletMolcas",processVerletMolcas), 
                     ("verletMolcasVel",processVerletMolcasVel),
                     ("molcasVel",processMolcasVel),("gaussVel",processGaussVel),
                     ("molcasTinker",processMolcasTinker),("molcasZeroVel",processMolcasZeroVelocity),                     
                     ("palmeiro",processPalmeiro), ("rewriteGateway",processGateway),
                     ("constrained",processConstrained),("prueba",processPrueba)]
                     
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

            
-- =========================>  Test API <=====================          
processPrueba :: Options -> IO ()
processPrueba opts =  do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[molcasFile] =  optInput opts      
  molcasInput  <- parseMolcasInputFile molcasFile
  print molcasInput
  
  
-- =============> Drivers to run the molecular dynamics simulations in Molcas <==============

processMolcas :: Options -> IO ()
processMolcas opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[xyz,molcasFile,input] =  optInput opts
  initData <- parseFileInput parseInput input
  let getter   = (initData ^.)      
  initialMol   <- initializeMolcasOntheFly xyz (getter getInitialState) temp
  molcasDriver getter opts molcasFile initialMol
  
processVerletMolcas :: Options -> IO ()
processVerletMolcas opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[xyz,molcasFile,input] =  optInput opts
  initData <- parseFileInput parseInput input
  let getter  = (initData ^.)
      project = getter getProject
  mol         <- initializeMolcasOntheFly xyz (getter getInitialState) temp
  molcasInput <- parseMolcasInputFile molcasFile
  let job     = Molcas molcasInput
  processVerlet getter opts job project mol

processVerletMolcasVel :: Options -> IO ()
processVerletMolcasVel opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[xyz,velxyz,molcasFile,input] =  optInput opts
  initData <- parseFileInput parseInput input
  let getter  = (initData ^.)
      project = getter getProject
  mol         <- initializeMolcasOntheFly xyz (getter getInitialState) temp
  vs         <- readInitialVel velxyz
  molcasInput <- parseMolcasInputFile molcasFile
  let job     = Molcas molcasInput
  processVerlet getter opts job project $ mol & getVel .~ vs   
    
processMolcasVel :: Options -> IO ()  
processMolcasVel opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[xyz,velxyz,molcasFile,input] =  optInput opts
  initData <- parseFileInput parseInput input
  let getter = (initData ^.)
  mol        <- initializeMolcasOntheFly xyz (getter getInitialState) temp
  vs         <- readInitialVel velxyz
  molcasDriver getter opts molcasFile $ mol & getVel .~ vs   

processMolcasZeroVelocity :: Options -> IO ()
processMolcasZeroVelocity opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[xyz,molcasFile,input] =  optInput opts
  initData <- parseFileInput parseInput input
  let getter   = (initData ^.)
      project  = getter getProject
  initialMol  <- initializeMolcasZeroVel xyz (getter getInitialState) temp
  molcasDriver getter opts molcasFile initialMol
  
molcasDriver :: (forall a. Getting a InitialDynamics a -> a) -> Options -> FilePath -> Molecule -> IO ()
molcasDriver getter opts molcasFile initialMol = do 
  let temp = fromMaybe 298 $ optTemperature opts
  molcasInput  <- parseMolcasInputFile molcasFile
  let numat         = initialMol ^. getAtoms . to length
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      thermo        = initializeThermo numat temp
      step          = 1
      job           = Molcas molcasInput
      aMatrix       = initialAMTX initialMol
      project  = getter getProject
  mol <- interactWith job project initialMol
  loggers <- mapM initLogger ["geometry.out", "result.out"]
  constantForceDynamics mol job thermo temp auTime audt (getter getForceAnchor) (getter getExtForceMod) aMatrix step project loggers
  mapM_ logStop loggers
        
-- | Molcas Tinker Interface         
processMolcasTinker :: Options -> IO ()
processMolcasTinker opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      files@[tinkerKey,tinkerXYZ,molcasFile,input] =  optInput opts      
  initData <- parseFileInput parseInput input
  let getter  = (initData ^.)
      project = getter getProject
  atomsQM               <- parserKeyFile tinkerKey
  molcasInput           <- parseMolcasInputFile molcasFile
  (initialMol,molcasQM) <- initializeMolcasTinker molcasFile (getter getInitialState) temp $ length atomsQM
  let numat         = initialMol ^. getAtoms . to length
      thermo        = initializeThermo numat temp
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      step          = 1
      job           = MolcasTinker molcasInput atomsQM molcasQM
  loggers <- mapM initLogger ["geometry.out", "result.out"]
  mol <- interactWith job project initialMol
  driverMolcasTinker mol audt auTime temp thermo job project step loggers
  mapM_ logStop loggers 

driverMolcasTinker :: Molecule -> Time -> DT -> Temperature -> Thermo -> Job -> Project -> Step -> [Logger] -> IO ()  
driverMolcasTinker mol t dt temp thermo job project step loggers = 
  if t <0 then return ()
          else do 
            let es = concatMap (printf "%.6f  ") $ mol ^. getEnergy . to head
            zipWithM_ ($) [printMol mol es, printData mol step] loggers
            (newMol,newThermo) <- dynamicNoseHoover mol dt temp thermo job project
            driverMolcasTinker newMol (t-dt) dt temp newThermo job project (succ step) loggers
            
-- =============> Functions to run the molecular dynamics simulations in Gaussian <==============

processVerletGaussian :: Options -> IO ()
processVerletGaussian opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      [input,fchk,out] = optInput opts 
  initData <- parseFileInput parseInput input
  let getter  = (initData ^.)
      project = "TullyExternalForces"
      theoryLevels  = getter getTheory
      basis         = getter getBasis
      job           = Gaussian (theoryLevels,basis)
  mol <- (updateMultiStates out) <=< (initializeSystemOnTheFly fchk $ getter getInitialState) $ temp    
  processVerlet getter opts job project mol

-- | on the fly molecular dynamics with applied external forces
processGauss :: Options -> IO ()
processGauss opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      [input,fchk,out] = optInput opts 
  initData <- parseFileInput parseInput input
  let getter = (initData ^.)
  mol <- (updateMultiStates out) <=< (initializeSystemOnTheFly fchk $ getter getInitialState) $ temp
  driverGaussian getter opts mol  

-- | user supplied initial velocities are used
processGaussVel :: Options -> IO ()
processGaussVel opts = do
  let temp = fromMaybe 298 $ optTemperature opts
      [input,fchk,out,velxyz] = optInput opts 
  initData   <- parseFileInput parseInput input
  let getter = (initData ^.)
  mol        <- (updateMultiStates out) <=< (initializeSystemOnTheFly fchk $ getter getInitialState) $ temp
  vs         <- readInitialVel velxyz
  driverGaussian getter opts $ mol & getVel .~ vs   
    
driverGaussian :: (forall a. Getting a InitialDynamics a -> a) -> Options -> Molecule -> IO ()
driverGaussian getter opts mol = do
  let temp = fromMaybe 298 $ optTemperature opts
      numat         = mol ^. getAtoms . to length
      thermo        = initializeThermo numat temp
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      aMatrix       = initialAMTX mol 
      step          = 1
      theoryLevels  = getter getTheory
      basis         = getter getBasis
      job           = Gaussian (theoryLevels,basis)
      project       = "TullyExternalForces"
  newMol  <- interactWith job project mol
  loggers <- mapM initLogger ["geometry.out", "result.out"] 
  constantForceDynamics newMol job thermo temp auTime audt (getter getForceAnchor) (getter getExtForceMod) aMatrix step project loggers
  mapM_ logStop loggers      
                                   
-- ==================> General Drivers <=================================  
  
processVerlet :: (forall a. Getting a InitialDynamics a -> a) -> Options -> Job ->  Project -> Molecule -> IO ()
processVerlet getter opts job project mol = do
  let numat         = mol ^. getAtoms . to length
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      aMatrix       = initialAMTX mol 
      step          = 1
  newMol  <- interactWith job project mol
  loggers <- mapM initLogger ["geometry.out", "result.out"] 
  driverVerlet newMol job auTime audt (getter getForceAnchor) (getter getExtForceMod) aMatrix step project loggers
  mapM_ logStop loggers

driverVerlet ::  Molecule -> Job -> Time -> DT-> Anchor -> Double -> MatrixCmplx -> Int -> Project  -> [Logger] -> IO () 
driverVerlet mol job time dt anchor externalForce aMatrix step project loggers = do
   if time < 0.0 then return ()
                 else do 
                  let es = concatMap (printf "%.6f  ") $ mol ^. getEnergy . to head
                  zipWithM_ ($) [printMol mol es, printData mol step] loggers                                         
                  newMol                <- velocityVerletForces mol dt job project anchor externalForce
                  (tullyMol,newAmatrix) <- tullyDriver dt aMatrix step newMol
                  printGnuplot newAmatrix tullyMol
                  let [oldRoot,newRoot] = (^.getElecSt) `fmap` [newMol,tullyMol]
                      newJob            =  if oldRoot == newRoot then job else updateNewJobInput job tullyMol
                  driverVerlet tullyMol newJob (time-dt) dt anchor externalForce newAmatrix (succ step) project loggers

  
constantForceDynamics ::  Molecule -> Job -> Thermo -> Temperature -> Time -> DT-> Anchor -> Double -> MatrixCmplx -> Int -> Project  -> [Logger] -> IO ()
constantForceDynamics mol job thermo temp time dt anchor externalForce aMatrix step project loggers = do
   if time < 0.0 then return ()
                 else do 
                  let es = concatMap (printf "%.6f  ") $ mol ^. getEnergy . to head
                  zipWithM_ ($) [printMol mol es, printData mol step] loggers                                         
                  (newMol,newThermo)    <- dynamicExternalForces mol dt temp thermo job project anchor externalForce
                  (tullyMol,newAmatrix) <- tullyDriver dt aMatrix step newMol
                  printGnuplot newAmatrix tullyMol
                  let [oldRoot,newRoot] = (^.getElecSt) `fmap` [newMol,tullyMol]
                      newJob            =  if oldRoot == newRoot then job else updateNewJobInput job tullyMol
                  constantForceDynamics tullyMol newJob newThermo temp (time-dt) dt anchor externalForce newAmatrix (succ step) project loggers

-- | Fewest Switches Tully Algorithm using Persico-Granucci Correction 
tullyDriver ::  DT -> MatrixCmplx -> Int ->  Molecule -> IO (Molecule,MatrixCmplx)
tullyDriver dt aMatrix step mol =
  if (mol^. getCoeffCI . to length) /= 3  -- At least 3 set of CI coefficients are required in oder to initialize the Tully
     then return (mol,aMatrix) 
     else tullyHS dt aMatrix step mol

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
  loggers <- mapM initLogger ["geometry.out", "result.out"]
  palmeiroLoop initialMol audt temp thermo job "" 1 loggers
  mapM_ logStop loggers       

  
palmeiroLoop :: Molecule -> DT -> Temperature -> Thermo -> Job -> String -> Step -> [Logger] -> IO ()
palmeiroLoop  mol dt t thermo job project step loggers = do
    print $ "Step: " ++ show step
    if t > 0 then do
                  (newMol,newThermo) <- dynamicNoseHoover mol dt t thermo job project  
                  zipWithM_ ($) [printMol mol "", printData mol step] loggers                                         
                  return ()
--                   palmeiroLoop newMol dt (t - dt) newThermo job project (succ step) loggers
             else return ()     
     
-- ===================> Miscallaneus Functions <================

-- | using a .key and .xyz tinker files a new Molcas input is generated
processGateway :: Options -> IO ()
processGateway opts =do
  let files@[tinkerKey,tinkerXYZ,molcasFile] =  optInput opts  
      (project,_)  = break (=='.') molcasFile
  atomsQM          <- parserKeyFile tinkerKey  
  tinkerQMMM       <- parserXYZFile tinkerXYZ
  molcasInput      <- parseMolcasInputFile molcasFile
  mol <- tinker2Molecule atomsQM tinkerQMMM defaultMol  
  let  numat       = length atomsQM
  molcasQM   <- parserInputMolcasQM molcasFile $ parserGatewayQM numat
  modifyMolcasInput molcasInput molcasQM project $ mol   
   
processConstrained :: Options -> IO ()
processConstrained = undefined


 

