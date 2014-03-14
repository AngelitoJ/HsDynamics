{-# Language DeriveFunctor ,BangPatterns, TupleSections #-}

module ClusterAPI (
                   launchCluster
                   ,writeShellPBS) where


import Control.Applicative
import Control.Concurrent (threadDelay)
import Control.Concurrent.Async
import Data.List as DL
import System.Directory 
import System.Exit(ExitCode( ExitFailure,ExitSuccess))
import System.Posix.User (getLoginName)
import System.Process (system,readProcess)
import Text.Parsec
import Text.Parsec.ByteString (parseFromFile)
import Text.Printf (printf)

-- Internal Imports
import CommonTypes

-- ==========================> <================================

-- | Command to summit jobs to a cluster  
launchCluster :: Command -> Args -> IO ()
launchCluster cmd name = do
  [pid] <- lines `fmap` readProcess cmd [name] ""
  simpleLoop pid 10 not    -- has the calculation started?
  simpleLoop pid 20 id -- has the calculation finished?
  
  
-- | check the cluster queue each t seconds
simpleLoop :: String  -> Int -> (Bool -> Bool) -> IO ()
simpleLoop pid t fun = do
  jobsId <- (concatMap words . lines) <$>  readProcess "qstat" ["-r"] ""
  if (fun $ DL.elem pid jobsId)
        then threadDelay (10^6*t) >> simpleLoop pid t fun -- Haskell use miliseconds
        else return ()

writeShellPBS :: Project ->  IO ()
writeShellPBS project = do
  hd   <- getCurrentDirectory
  name <- getLoginName
  let   file    = project ++ ".sh"        
        x1  =  "#!/bin/bash -l"
        x2  =  "#PBS -N " ++ project
        x3  =  "#PBS -S /bin/sh"
        x4  =  "#PBS -r n"
        x5  =  "#PBS -l nodes=1:ppn=4"
        x6  =  "#PBS -l mem=8000MB"
        x7  =  "rand=$(basename $PBS_JOBID .master)"
        x8  =  "source /apps/intel/composer_xe_2011_sp1.9.293/bin/compilervars.sh intel64"
        x9  =  "export Project=" ++ project
        x10 =  "export MOLCAS=/apps/Chem/molcas7.9"
        x11 =  "export MOLCASMEM=8000MB"
        x12 =  "export TINKER=/apps/Chem/molcas7.9/tinker/bin/" 
        x13 =  "HomeDir=" ++ hd 
        x14 =  "cd /scratch/" ++ name 
        x15 =  "mkdir " ++ project ++ ".$rand" 
        x16 =  "TempDir=/scratch/" ++ name ++ "/" ++ project ++ ".$rand/" 
        x17 =  "export WorkDir=$TempDir"
        x18 =  "export InpDir=$HomeDir"
        x19 =  "cp $HomeDir/* $WorkDir/"
        x20 =  "cd $WorkDir"
        x21 =  "ln -sf /apps/bin/molcas molcas"
        x22 =   "./molcas $WorkDir/$Project.input >$InpDir/$Project.out 2>$InpDir/$Project.err"
        x23 =   "wait"        
  writeFile file $ concatMap (++"\n") [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23] 
  p <- getPermissions file
  setPermissions file (p {executable = True})

    
  


