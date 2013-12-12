
module TinkerQMMM (
                   parserKeyFile
                  ,parserXYZFile
                  ,reWriteXYZtinker
                   )  where

import Control.Applicative
import Control.Lens
import qualified Data.Array.Repa as R
import Data.List.Split (chunksOf)
import qualified Data.Vector as V
import System.Directory (renameFile)
import Text.Parsec
import Text.Parsec.ByteString
import Text.Printf

-- =========> Internal imports <========
import CommonTypes 
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
                   return (l,n)
  
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

reWriteXYZtinker :: Molecule -> [(Label,Int)] -> Command -> Project ->  IO ()
reWriteXYZtinker mol atomsQM commandTinker project = do
   renameFile (project++".xyz") "temp"
   tinkerQMMM <- parserXYZFile "temp"
   let coordinates = chunksOf 3 $ R.toList $ mol^.getCoord
       numbersQM   = snd <$> atomsQM
       dim         = pred . length $ numbersQM
       qmmm        = updateAtomsMM numbersQM 
       updateAtomsMM xs = fst $ foldr step ([],reverse coordinates) tinkerQMMM      
       step atom t@(acc,(x:xs)) =                      -- if the atom is QM append the new Coordinates
            if  (atom^.numberMM) `elem` numbersQM      -- else return the same atom
                  then let newAtom = set xyzMM x atom  -- the initial acc has the coordinates in reversed
                       in (newAtom:acc,xs)             -- order because we are replacing the QM atoms beginning
                  else t                               -- from the last one
   writeXYZtinker qmmm  project  
 
writeXYZtinker :: [AtomMM] -> Project -> IO ()
writeXYZtinker atomsMM project = do
  let numat = printf "  %4d\n" $ length atomsMM
      dat = concatMap (\(AtomMM number label [x,y,z] args) -> printf "%5d %s %11.6f %11.6f %11.6f %s\n" number label x y z args) atomsMM
  writeFile (project ++ ".xyz") $ numat ++ dat