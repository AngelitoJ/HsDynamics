module Paths_HsDynamics (
    version,
    getBinDir, getLibDir, getDataDir, getLibexecDir,
    getDataFileName
  ) where

import qualified Control.Exception as Exception
import Data.Version (Version(..))
import System.Environment (getEnv)
catchIO :: IO a -> (Exception.IOException -> IO a) -> IO a
catchIO = Exception.catch


version :: Version
version = Version {versionBranch = [0,1,0,0], versionTags = []}
bindir, libdir, datadir, libexecdir :: FilePath

bindir     = "/home/felipe/.cabal/bin"
libdir     = "/home/felipe/.cabal/lib/HsDynamics-0.1.0.0/ghc-7.4.2"
datadir    = "/home/felipe/.cabal/share/HsDynamics-0.1.0.0"
libexecdir = "/home/felipe/.cabal/libexec"

getBinDir, getLibDir, getDataDir, getLibexecDir :: IO FilePath
getBinDir = catchIO (getEnv "HsDynamics_bindir") (\_ -> return bindir)
getLibDir = catchIO (getEnv "HsDynamics_libdir") (\_ -> return libdir)
getDataDir = catchIO (getEnv "HsDynamics_datadir") (\_ -> return datadir)
getLibexecDir = catchIO (getEnv "HsDynamics_libexecdir") (\_ -> return libexecdir)

getDataFileName :: FilePath -> IO FilePath
getDataFileName name = do
  dir <- getDataDir
  return (dir ++ "/" ++ name)
