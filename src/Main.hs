{-# Language RankNTypes #-}
-- HsDynamics: 

-- @2013 Felipe Zapata, Alessio Valentini, Angel Alvarez from The ResMol Group  

module Main where


import Data.List (transpose,zipWith3)
import Data.Complex
import Data.Maybe ( fromMaybe )
import qualified Data.Vector.Unboxed as VU
import Control.Applicative
import Control.Concurrent
import Control.Concurrent.Async
import Control.Parallel.Strategies (parMap,rdeepseq)
import Control.Lens ((.~),(^.),(&),Getting(..),to)
import Control.Monad ((<=<),liftM,zipWithM,zipWithM_)
import Control.Monad.Trans.Either
import System.Directory (createDirectoryIfMissing)
import System.Environment ( getArgs )
import System.FilePath
import System.Cmd ( system )
import System.Console.GetOpt
import Text.Printf

import Control.Exception (Exception(..),catch)
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
                     ("verletMolcasVel",processVerletMolcasVel),("verletGaussVel",processVerletGaussVel),
                     ("NVTMolcas",processNVTMolcas),
                     ("molcasZeroVel",processMolcasZeroVelocity),                                          
                     ("molcasVel",processMolcasVel),("gaussVel",processGaussVel),
                     ("molcasTinker",processMolcasTinker), ("molcasTinkerVel",processMolcasTinkerVel),                     
                     ("palmeiro",processPalmeiro), ("rewriteGateway",processGateway),
                     ("readOut",processReadOut),("calcInternal",processCalcInternals),
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
processPrueba opts = undefined

-- ======> Share Data and Functions <=================

topData :: Options -> IO TopData
topData opts =  do
  let temp      = fromMaybe 298 $ optTemperature opts
      (input:_) =  optInput opts
  initData   <- parseFileInput parseInput input
  let getter  = (initData ^.)
  return $ TopData  getter temp
  
-- =============> Drivers to run the molecular dynamics simulations in Molcas <==============

processMolcas :: Options -> IO ()
processMolcas opts = molcas opts Nothing
  
processMolcasVel :: Options -> IO ()  
processMolcasVel opts = do
  let [input,xyz,velxyz,molcasFile] =  optInput opts
  vs  <- readInitialVel velxyz
  molcas opts $ Just vs   
  
processMolcasZeroVelocity :: Options -> IO ()
processMolcasZeroVelocity opts = do
  let [input,xyz,molcasFile] =  optInput opts
  vs <- initializeMolcasZeroVel xyz
  molcas opts $ Just vs

-- | Molcas Tinker Interface           
processMolcasTinker :: Options -> IO ()
processMolcasTinker opts = molcasTinker opts Nothing

processMolcasTinkerVel :: Options -> IO ()
processMolcasTinkerVel opts = do
  let [input,tinkerKey,molcasFile,velxyz] =  optInput opts      
  vs  <- readInitialVel velxyz    
  molcasTinker opts $ Just vs     
  
molcas :: Options -> Maybe Vec -> IO ()
molcas opts vel =  do
  let (input:xyz:molcasFile:_) =  optInput opts
  TopData getter temp <-  topData opts
  molcasInput         <- parseMolcasInputFile molcasFile
  let job             = Molcas molcasInput
  initialMol          <- initializeMolcasOntheFly xyz (getter getInitialState) temp
  molcasDriver getter opts job temp $ addExternalVec initialMol vel
  
molcasTinker :: Options -> Maybe Vec -> IO ()
molcasTinker opts vel = do
  let (input:tinkerKey:molcasFile:_) =  optInput opts      
  TopData getter temp   <- topData opts    
  atomsQM               <- parserKeyFile tinkerKey
  createDirectoryIfMissing False "TinkerGeometries" -- create directory to save the .xyz generated by tinker
  molcasInput           <- parseMolcasInputFile molcasFile
  (initialMol,molcasQM) <- initializeMolcasTinker molcasFile (getter getInitialState) temp $ length atomsQM
  molcasDriver getter opts (MolcasTinker molcasInput molcasQM) temp $ addExternalVec initialMol vel

molcasDriver :: (forall a. Getting a InitialDynamics a -> a) -> Options -> Job -> Temperature -> Molecule -> IO ()
molcasDriver getter opts job temp initialMol = do 
  let numat         = initialMol ^. getAtoms . to length
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      thermo        = initializeThermo numat temp
      step          = 1     
      aMatrix       = initialAMTX initialMol
      project  = getter getProject
  mol <- firsStepMolcas job project step initialMol
  loggers <- mapM initLogger ["geometry.out", "result.out","totalEnergy.out"]
  constantForceDynamics mol job thermo temp auTime audt (getter getForceAnchor) (getter getExtForceMod) aMatrix step project loggers
  mapM_ logStop loggers
    
-- | Simple Constant temperature dynamics
processNVTMolcas :: Options -> IO ()
processNVTMolcas opts =do
  let [input,xyz,molcasFile] =  optInput opts
  TopData getter temp   <- topData opts          
  molcasInput  <- parseMolcasInputFile molcasFile
  let project  = getter getProject
      job      =  Molcas molcasInput
  initialMol   <- initializeMolcasOntheFly xyz (getter getInitialState) temp
  driverNVT getter opts job project temp initialMol
             
-- | Velocity Verlet Molcas                           
processVerletMolcas :: Options -> IO ()
processVerletMolcas opts = verletMolcas opts Nothing
  
processVerletMolcasVel :: Options -> IO ()
processVerletMolcasVel opts = do
  let [input,xyz,molcasFile,velxyz] =  optInput opts
  vs   <- readInitialVel velxyz
  verletMolcas opts $ Just vs   

verletMolcas :: Options -> Maybe Vec -> IO ()
verletMolcas opts vel =  do
  let (input:xyz:molcasFile:_) =  optInput opts
  TopData getter temp   <- topData opts                
  let project = getter getProject
  mol         <- initializeMolcasOntheFly xyz (getter getInitialState) temp
  molcasInput <- parseMolcasInputFile molcasFile
  let job     = Molcas molcasInput
  processVerlet getter opts job project $ addExternalVec mol vel    
                               
-- =============> Functions to run the molecular dynamics simulations in Gaussian <==============
processGauss :: Options -> IO ()
processGauss opts = gaussian opts Nothing

-- | on the fly molecular dynamics with applied external forces
gaussian :: Options -> Maybe Vec -> IO ()
gaussian opts vel = do
  let (input:fchk:out:_) = optInput opts 
  TopData getter temp   <- topData opts          
  mol <- (updateMultiStates out fchk) <=< (initializeSystemOnTheFly fchk $ getter getInitialState) $ temp
  driverGaussian getter opts $ addExternalVec mol vel

processVerletGaussian :: Options -> IO ()
processVerletGaussian opts = verletGaussian opts Nothing

processVerletGaussVel :: Options -> IO ()
processVerletGaussVel opts = do
  let [input,fchk,out,velxyz] =  optInput opts
  vs  <- readInitialVel velxyz
  verletGaussian opts $ Just vs

verletGaussian :: Options -> Maybe Vec -> IO () 
verletGaussian opts vel = do
  let (input:fchk:out:_) = optInput opts 
  TopData getter temp   <- topData opts          
  let project           = "TullyExternalForces"
      theoryLevels      = getter getTheory
      basis             = getter getBasis
      job               = Gaussian (theoryLevels,basis)
  mol <- (updateMultiStates out fchk) <=< (initializeSystemOnTheFly fchk $ getter getInitialState) $ temp    
  processVerlet getter opts job project $ addExternalVec mol vel

-- | user supplied initial velocities are used
processGaussVel :: Options -> IO ()
processGaussVel opts = do
  let [input,fchk,out,velxyz] = optInput opts 
  vs  <- readInitialVel velxyz
  gaussian opts $ Just vs    
  
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
  newMol  <- interactWith job project step mol
  loggers <- mapM initLogger ["geometry.out", "result.out","totalEnergy.out"]
  constantForceDynamics newMol job thermo temp auTime audt (getter getForceAnchor) (getter getExtForceMod) aMatrix step project loggers
  mapM_ logStop loggers      
      
-- ==================> General Drivers <=================================  

driverVerletGround :: (forall a. Getting a InitialDynamics a -> a) -> Options -> Job ->  Project -> Molecule -> IO ()
driverVerletGround getter opts job project mol = do
  let numat         = mol ^. getAtoms . to length
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      step          = 1
  loggers <- mapM initLogger ["geometry.out", "result.out","totalEnergy.out"]
  loop loggers job auTime audt step project mol
              
  where
    loop logs job time dt step project molecule = 
         if time < 0 then mapM_ logStop logs
                     else do
                        newMol  <- velocityVerletForces molecule dt job project [] 0 step
                        let es = concatMap (printf "%12.6f  ")  $ newMol ^. getEnergy . to head
                        zipWithM_ ($) [printMol newMol es, printData newMol step, printTotalEnergy newMol] logs
                        loop logs job (time-dt) dt (succ step) project newMol
  
-- | Loop to run Molecular dynamics at constant temperature without using the Tully algorithm  neither external forces
driverNVT :: (forall a. Getting a InitialDynamics a -> a) -> Options -> Job ->  Project -> Temperature -> Molecule -> IO ()
driverNVT getter opts job project temp mol = do
  let numat         = mol ^. getAtoms . to length
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      step          = 1
      bath     = initializeThermo numat temp      
  loggers <- mapM initLogger ["geometry.out", "result.out","totalEnergy.out"]
  loop loggers job auTime audt step project temp bath mol
  
  where
    loop logs job time dt step project temp bath molecule = 
         if time < 0 then mapM_ logStop logs
                     else do
                        (newMol,newThermo) <- dynamicNoseHoover molecule dt temp bath job project step 
                        let es = concatMap (printf "%12.6f  ")  $ newMol ^. getEnergy . to head
                        zipWithM_ ($) [printMol newMol es, printData newMol step, printTotalEnergy newMol] logs
                        loop logs job (time-dt) dt (succ step) project temp newThermo newMol  
  

processVerlet :: (forall a. Getting a InitialDynamics a -> a) -> Options -> Job ->  Project -> Molecule -> IO ()
processVerlet getter opts job project mol = do
  let numat         = mol ^. getAtoms . to length
      [auTime,audt] = fmap (/au_time) $ getter `fmap` [getTime,getdt] 
      aMatrix       = initialAMTX mol 
      step          = 1
  newMol  <- interactWith job project step mol
  loggers <- mapM initLogger ["geometry.out", "result.out","totalEnergy.out"] 
  driverVerlet newMol job auTime audt (getter getForceAnchor) (getter getExtForceMod) aMatrix step project loggers
  mapM_ logStop loggers

driverVerlet ::  Molecule -> Job -> Time -> DT-> Anchor -> Double -> MatrixCmplx -> Int -> Project  -> [Logger] -> IO () 
driverVerlet mol job time dt anchor externalForce aMatrix step project loggers = do
   if time < 0 then return ()
               else do 
                let es = concatMap (printf "%.6f  ") $ mol ^. getEnergy . to head
                zipWithM_ ($) [printMol mol es, printData mol step, printTotalEnergy mol] loggers                                         
                newMol                <- velocityVerletForces mol dt job project anchor externalForce step
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
                  zipWithM_ ($) [printMol mol es, printData mol step, printTotalEnergy mol] loggers
                  (newMol,newThermo)    <- dynamicExternalForces mol dt temp thermo job project anchor externalForce step
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
      [input,xyz] =  optInput opts         
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
                  (newMol,newThermo) <- dynamicNoseHoover mol dt t thermo job project step 
                  zipWithM_ ($) [printMol mol "", printData mol step] loggers                                         
                  return ()
--                   palmeiroLoop newMol dt (t - dt) newThermo job project (succ step) loggers
             else return ()     
     
-- ===================> Miscallaneus Functions <================
addExternalVec :: Molecule -> Maybe Vec -> Molecule
addExternalVec mol vel = 
  case vel  of
       Nothing -> mol
       Just vs -> mol & getVel .~ vs

-- | using a .key and .xyz tinker files a new Molcas input is generated
processGateway :: Options -> IO ()
processGateway opts =do
  let [tinkerKey,tinkerXYZ,molcasFile] =  optInput opts  
      (project,_)  = break (=='.') molcasFile
  atomsQM          <- parserKeyFile tinkerKey  
  tinkerQMMM       <- parserXYZFile tinkerXYZ
  molcasInput      <- parseMolcasInputFile molcasFile
  mol <- tinker2Molecule atomsQM tinkerQMMM defaultMol  
  let  numat       = length atomsQM
  molcasQM   <- parserInputMolcasQM molcasFile $ parserGatewayQM numat
  modifyMolcasInput molcasInput molcasQM project $ mol   

processReadOut :: Options -> IO ()
processReadOut opts = do
  let file@[out] =  optInput opts 
      state          = S1
  initialMol  <- initializeMolcasOntheFly out state 300
  mols <- parserGeomVel out initialMol         
      {-totalS0 = parMap rdeepseq calcTotalEnergy $ (& getElecSt .~ Left S0) `fmap` mols-}
      {-total   = parMap rdeepseq calcTotalEnergy mols-}
  let fun  (x:xs)       = (x,if null xs then 0 else head xs)
      (totalS0,totalS1) = unzip $ parMap rdeepseq fun $ (^. getEnergy . to head) `fmap` mols
      kinetic           = parMap rdeepseq (\x -> calcEk (x ^. getVel) (x ^. getMass)) mols
      a1                = writeFile ("EnergyS0") $ concatMap (printf "%.5f\n") totalS0
      a2                = writeFile ("EnergyS1") $ concatMap (printf "%.5f\n") totalS1
      a3                = writeFile ("KineticEnergy") $ concatMap (printf "%.5f\n") kinetic
  parallelLaunch [a1,a2,a3]
  print "Done"

processCalcInternals :: Options -> IO ()
processCalcInternals opts = do
  let file@[fileConex,out] =  optInput opts
      state                =  S1
  conex      <- parserFileInternas fileConex
  initialMol <- initializeMolcasOntheFly out state 300
  mols       <- parserGeomVel out initialMol
  let numat           = initialMol ^. getAtoms . to length
      numberBonds     = pred numat
      internals       = parMap rdeepseq (calcInternals conex) mols
      funSlice        = VU.toList . VU.slice 0 numberBonds 
      funPrint        = concatMap (printf "%12.6f\n")
      fun             = transpose 
      funWrite i ints = async . writeFile ("enlace_" ++ (show i) ++ ".out") . funPrint $ ints
  ids <- zipWithM (funWrite) [1..] $ transpose . fmap funSlice $ internals
  mapM_ wait ids
  
processConstrained :: Options -> IO ()
processConstrained = undefined


 