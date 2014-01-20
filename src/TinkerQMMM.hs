{-# Language BangPatterns, TupleSections #-}

module TinkerQMMM (
                   parserKeyFile
                  ,parserXYZFile
                  ,reWriteXYZtinker
                  ,tinker2Molecule
                   )  where

import Control.Applicative
import Control.Exception (bracket)
import Control.Lens
import Data.Array.Repa as R hiding ((++))
import Data.List (lookup,unfoldr)
import Data.List.Split (chunksOf)
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
import System.Directory (renameFile,removeFile)
import Text.Parsec
import Text.Parsec.ByteString
import Text.Printf

-- =========> Internal imports <========
import CommonTypes 
import Constants (a0)
import Molcas
import ParsecNumbers
import ParsecText


-- ===========> Types <===========

-- ===========> Parse .key tinker <===========
parserKeyFile :: FilePath -> IO [(Label,Int)]
parserKeyFile filename = do 
  r <- parseFromFile parseKeyTinker filename
  case r of 
     Left msg -> error $ show msg
     Right xs -> return xs

-- REMEMBER that
-- data ParsecT s u m a
-- type Parsec s u = ParsecT s u Identity
-- Parser = Parsec ByteString ()

parseKeyTinker ::  Parser [(Label,Int)]
parseKeyTinker = do 
               anyLine --header
               (_,numat) <- parserLabelNumat -- number of QM atoms
               count numat parserLabelNumat
  
parserLabelNumat :: Parser (Label,Int)
parserLabelNumat = do 
                   l <- manyTill anyChar space
                   spaces 
                   n <- intNumber
                   anyLine
                   return (l, n)
  
-- =========================> Parser Tinker xyz <===============

-- | Parse the XYZ in the Tinker format 
parserXYZFile :: FilePath -> IO [AtomMM]
parserXYZFile filename = do
  r <- parseFromFile parseXYZ filename
  case r of 
       Left msg -> error $ show msg
       Right xs -> return xs

-- | First count the number of atoms the parse line by line
parseXYZ :: Parser [AtomMM]
parseXYZ = do
           numat <- (read . head. words) <$> anyLine
           count numat parseLineAtomMM 
                      
parseLineAtomMM :: Parser AtomMM
parseLineAtomMM = do 
                 numat <- (spaces >> intNumber)
                 label <- (spaces >> manyTill anyChar space)
                 manyTill anyChar space 
                 xs    <- count 3 (spaces >> realNumber)
                 ps    <- anyLine
                 return $ AtomMM numat label xs ps
                 
-- ======================> rewrite XYZ File <============

reWriteXYZtinker :: Molecule -> [(Label,Int)] -> Project ->  IO ()
reWriteXYZtinker mol atomsQM project = 
   bracket (renameFile (project++".xyz") "temp") (\_ -> removeFile "temp") $ \_ -> do
   tinkerQMMM <- parserXYZFile "temp"
   let coordinates   = chunksOf 3 $ R.toList . R.computeUnboxedS . R.map (*a0) $ mol^.getCoord -- Coordinates are printed in Amstrongs
       numbersQM     =  snd <$> atomsQM  
       qmmm          = updateAtomsQM
       mapa          = zip numbersQM coordinates 
       fun           = fromMaybe (error "unknown atomic numbersQM") . flip lookup mapa
       updateAtomsQM = fmap step tinkerQMMM                -- if the atom is QM append the new Coordinates
       step atom     = let indx = atom^.numberMM           -- else return the same atom
                           newPosition = fun indx
                       in if indx `elem` numbersQM
                          then atom &  xyzMM .~ newPosition
                          else atom       
   writeXYZtinker qmmm project          
              
writeXYZtinker :: [AtomMM] -> Project -> IO ()
writeXYZtinker atomsMM project = do
  let numat = printf "  %4d\n" $ length atomsMM
      dat = concatMap (\(AtomMM number label [x,y,z] args) -> printf "%5d %s %11.6f %11.6f %11.6f %s\n" number label x y z args) atomsMM
  writeFile (project ++ ".xyz") $ numat ++ dat
  
-- ============> <===============

tinker2Molecule :: [(Label,Int)] -> [AtomMM] -> Molecule -> IO Molecule
tinker2Molecule atomsQM  tinkerQMMM  mol= do
  let  numat       = length atomsQM
       numbersQM   = (snd) <$> atomsQM -- tinker indexes begin at 1 therefore indexes are traslating using pred
       lookQM atom = (atom^.numberMM) `elem` numbersQM 
       coords      = R.fromListUnboxed (Z:. 3*numat) . concatMap (^.xyzMM) $ filter lookQM tinkerQMMM  
  return $ mol & getCoord .~ coords
  