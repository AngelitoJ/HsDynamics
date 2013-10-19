
-- Monadic command line checkers..
-- @2013 Angel Alvarez

module OptsCheck where

-- Main module needed imports
import Control.Exception (SomeException,evaluate,try)
import Control.Monad(foldM,liftM,ap)
import Control.Monad.IO.Class
import Control.Monad.Trans.Either
import Data.List (find)
import Data.Maybe ( fromMaybe )
import Data.Either
import Text.Show.Functions
import System.Environment
import System.Console.GetOpt
import System.Directory ( doesDirectoryExist, doesFileExist )
import System.FilePath

-- Record for storing cmdline options
data Options = Options
 { optDump        :: Bool
 , optModules     :: [(String,(Options-> IO()))]
 , optMode        :: Maybe (Options->IO())
 , optVerbose     :: Bool
 , optShowVersion :: Bool
 , optOutput      :: Maybe FilePath
 , optDataDir     :: Maybe FilePath
 , optInput       :: [FilePath]
 , optTemperature :: Maybe Double
 } deriving (Show)

        
-- An EitherT container to store parsed opts from commandline or error messages
type OptsResult = EitherT String IO Options

-- An Opts filter runnng in the EitherT IO Stack
type OptsFilter = ( Options -> OptsResult )

-- A Policy describing a command line options with a checking filter
type OptsPolicy = OptDescr OptsFilter


-- =============================================== Options checking ==============================================
-- progOpts: getOpt args processor that also checks semantically getopt results or bang in front of the user
-- Upon checking with getopt this function gets a list of lambdas representing semantical checks
-- on every cmdline switch, as an example; we check for input file presence and the check that either data
-- director is a valid one and input file exists and is readable. Those checks are performed by
-- filename_check, check_datadir and check_input_file respectively. These funs are stored
-- in the Options structure that getopt uses.

--  We use a EitherT transformer to combine Either chaining with arbitrary IO actions needed during checking
-- ===============================================================================================================

progOpts :: [String] -> Options -> [OptsPolicy] -> OptsResult
progOpts args defaultOptions acceptedOptions =
   case getOpt RequireOrder acceptedOptions args of
      (funs,[],[]) -> do
          left "input file(s) missing"
      (funs,filenames,[]) -> do
          resultOfFuns <- foldl (>>=) (return defaultOptions) funs               -- Perform monadic checkings upon getOpt supplied functions
          foldM check_input_file resultOfFuns $ reverse filenames                          -- Now check if all the input files exist and are accesible
      (_,_,errs) -> do
          left ( concat errs )


-- =============================================== Monadic Options checkers =======================================
-- getOpt will partially apply against the supplied argument do we can just come over the options record
          
          
-- Who we are?, sort of alien outbreak?
check_version :: Options -> OptsResult
check_version optsR = return $ optsR { optShowVersion = True }
          
-- help message, dont panic, we are just not going anywhere after showing the help
check_help :: Options -> OptsResult
check_help _ = left "Command line help requested"
    
--check supplied input files exist or bang if any is not.
check_input_file :: Options -> String -> OptsResult
check_input_file optsR@Options { optInput = files , optDataDir = dataDir } filename = do
    test <- liftIO $ filename_check dataDir filename
    case test of
         True -> return $ (optsR { optInput = filename : files })
         False -> left $ "input file "++ filename ++ " not readable"
    where
        filename_check :: Maybe FilePath -> FilePath -> IO Bool --check file with or without data directory
        filename_check (Just datadir) filename = doesFileExist $ combine datadir filename
        filename_check Nothing filename        = doesFileExist filename

-- User passed some directory, we make sure this dir exits so file will matched against it
check_data_dir :: String -> Options -> OptsResult
check_data_dir dir optsR = do
    test <- liftIO $ doesDirectoryExist dir
    case test of
         True -> return $ optsR { optDataDir =  Just dir } 
         False -> left ( "Data directory " ++ dir ++ " does not exist" )

-- check user wants verbosity level increased
check_verbosity :: Options -> OptsResult
check_verbosity optsR = return $ optsR { optVerbose = True }

-- check mode of operation. A list of modules is provided in the options record
check_operation_mode :: String -> Options -> OptsResult 
check_operation_mode mode optsR@Options { optModules = modules } = do
    return $ optsR { optMode = selectedModule }
    where
        selectedModule = case (findmodule mode modules) of
                              Just (_,fun) -> Just fun
                              Nothing      -> Nothing
        findmodule :: String -> [(String, (Options-> IO()))] -> Maybe (String,(Options -> IO ()))
        findmodule mode = find ((== mode ).fst )


-- dump either options or errors as we get passthrought
check_dump_options :: Options -> OptsResult
check_dump_options optsR = do
    liftIO $ putStrLn $ "\n\nOptions dumping selected record: \n\t" ++ show optsR ++ "\n"
    return $ optsR { optDump =True } 
     
-- | Optional temperature argument take from the command line
check_temperature :: String -> Options -> OptsResult
check_temperature temp optsR = do
    r <- liftIO (try (evaluate (read temp)) :: IO (Either SomeException Double))
    case r of
         Right val -> if val < 0 then error "Temperature must be a positive number"
                                 else return $ optsR { optTemperature = Just val }
         Left _ -> left $ "Wrong temperature value: " ++ temp   

