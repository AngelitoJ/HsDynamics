-- ------------------------------------------------------------------------------------
--  Gaussian 2003 checkpointing file parser
--  Authors: Angel Alvarez, Felipe Zapata
--
-- data file consists on two lined header and one or more data blocks
-- every item in a line is separated by one or more spaces
-- There are single line datablocks 'S' or a multi line datablocks 'M'
-- All data block start with a one line header 
-- 'S' blocks contain a text label, a type field and a single value
-- 'M' blocks contain a text label, a type field and cardinality field
-- A data type filed consist one of 'I' , 'R' for for Integer data , real data etc...
-- A cardinality field consists of 'N='  and a integer representing amount of data
-- Integer values appear in lines, up to six values per line
-- Real values appear in lines up to five values per line

-- Some header types...
-- "Number of atoms                            I               26"
-- "Info1-9                                    I   N=           9"
-- "Nuclear charges                            R   N=          26"
-- "Int Atom Types                             I   N=          26"
-- ------------------------------------------------------------------------------------

module Gaussian where

import Data.Char
import Data.Function (on)
import Data.List (foldl',sortBy)
import Text.ParserCombinators.Parsec
import Text.Parsec.Numbers

import CommonTypes 

-- ===============> Show Instances <================

showRblock :: String -> Int -> [Double] -> String
showRblock l i v 
    | i == 1 = "Real block : " ++ (show l) ++ " Element : " ++ (show v)
    | i <= 5 = "Real block : " ++ (show l) ++ " Elements : " ++ (show v)
    | otherwise = "Real block : " ++ (show l) ++ "\n\tElements(" ++ (show i) ++ "): " ++ (show v)

showIblock :: String -> Int -> [Integer] -> String
showIblock  l i v 
    | i == 1 = "Integer block : " ++ (show l) ++ " Element : " ++ (show v)
    | i <= 6 = "Integer block : " ++ (show l) ++ " Elements : " ++ (show v)
    | otherwise = "Integer block : " ++ (show l) ++ "\n\tElements(" ++ (show i) ++ "): " ++ (show v)

showTblock :: String -> String
showTblock l = "Text    block : " ++ (show l)

instance Show GauBlock where
    show (IGauBlock l i v) = showIblock l i v
    show (RGauBlock l i v) = showRblock l i v
    show (TGauBlock l)     = showTblock l 


-- Predicate Functions 
isGaussGrad (RGauBlock label n xs) = let labelGrad = ["Cartesian","Gradient"]
                                         p = (words label) == labelGrad
                                     in  if p then True else False 


--  ================> FCHK Parser <=====================

-- Note :
-- "Return" refers to lifting results in the proper monad
-- as most consumers are monadic, most producers just "return" 
-- data in the general sense

-- | Parse the type of gaussian computation
readTheoryLevel :: String -> IO String
readTheoryLevel file = do
  r <- parseFromFile parseTheoryLevel file
  case r of
       Left msg -> error . show $ msg
       Right xs -> return xs
       
parseTheoryLevel :: GenParser Char st String
parseTheoryLevel = do
  manyTill anyChar $ char '#'
  r <- manyTill anyChar $ (try $ string "Charge")
  return  $ "#" ++ r
  
-- Parse a file and returns either an error string or a list of data blocks

parseGaussianCheckpoint :: String -> IO (Either ParseError [GauBlock])
parseGaussianCheckpoint = parseFromFile gaussianChkParser

-- Parse Gaussian Checkpoint format
gaussianChkParser :: GenParser Char st [GauBlock]
gaussianChkParser = do
    fileheader <- count 2 anyline      -- discard first two lines of header
    datablocks <- many1 goodDataBlock  -- parse at least one datablock
    eof
    return datablocks              -- 

-- Parse a text line as a whole into a string
anyline :: GenParser Char st String
anyline = manyTill anyChar newline     -- whatever we found till we hit a newline

--sanitize datablocks
goodDataBlock :: GenParser Char st GauBlock
goodDataBlock = do
    candidate <- datablock
    isGood <- blockcheck candidate   -- check for sane blocks
    if isGood
       then return candidate
       else fail "GauBlock cardinality is not good"  -- TODO: Just tell us what was the offeding block

-- Try to parse one of the known datablock types or give up
datablock :: GenParser Char st GauBlock
datablock = do
    try singleblock    -- Is it a single block?
    <|> multiblock     -- can it be multi block?
    <|> whatever       -- Well, sure we can put it in a text block
    <?> "Cachis!"      -- Something gone wrong here we cannot not parse anything at all!

   -- Parse while text lines into TGauBlocks   
whatever :: GenParser Char st GauBlock
whatever = do
    iFound <- anyline
    return $ TGauBlock iFound

-- Check block have proper number of elements
blockcheck :: GauBlock -> GenParser Char st Bool
blockcheck (IGauBlock t n v) = return $ n == length v
blockcheck (RGauBlock t n v) = return $ n == length v
blockcheck (TGauBlock t) = return True


-- Try to parse single data block that contains a header and one datavalue somewhat spaced
singleblock :: GenParser Char st GauBlock
singleblock = try $ do
    (blockLabel,valueType) <- simpleheader
    case valueType of
        'I' -> do
                skipMany1 space
                valueData <- parseIntegral
                newline
                return $ IGauBlock blockLabel 1 [valueData]
        'R' -> do
                skipMany1 space
                valueData <- parseFloat
                newline
                return $ RGauBlock blockLabel 1 [valueData]

-- Try to parse a multi block containing a header and two or more values across several lines (up to five or six elements every line)
multiblock :: GenParser Char st GauBlock
multiblock = try $ do
    (blockLabel,valueType,valueNumber) <- multiheader
    case valueType of
        'I' -> do
                values              <- ivalues valueNumber
                return $ IGauBlock blockLabel valueNumber values           -- IGauBlocks store integer data
        'R' -> do
                values              <- rvalues valueNumber
                return $ RGauBlock blockLabel valueNumber values           -- RGauBlocks store real data

-- Parse simple block header consisting of a label, and a type
simpleheader :: GenParser Char st (String,Char)
simpleheader = do
    dlabel <- labelfield
    dtype  <- oneOf "IR"         -- data types are 'I'nteger or 'R'eal
    return $ (dlabel,dtype)

-- Parse a multi block header line consisting of a label, a type, and a cardinality
multiheader :: GenParser Char st (String,Char,Int)
multiheader = do
    dlabel <- labelfield
    dtype  <- oneOf "IR"         -- data types are 'I'nteger or 'R'eal
    spaces
    dcar   <- cardinality        -- How many values will have to parse?
    newline
    return $ (dlabel,dtype,dcar)

-- Parse the label (always 43 chars)
labelfield :: GenParser Char st String
labelfield = do
    result <- count 43 anyChar
    return $ result

-- ====================> Log Parser <===================================
-- Parse the Information which is not contain in the fchk file

  
-- Parse the Energies and coeffcient of CASSCF calculations
  
readEigenDataCASSCF :: FilePath -> IO (Either ParseError [EigenBLock])
readEigenDataCASSCF = parseFromFile parseEigenData

parseLogGaussian :: Int -> FilePath -> IO (Either ParseError GaussLog)
parseLogGaussian k = parseFromFile combination
  where combination = do
        eg   <- parseEigenData
        grad <- parseLogGrad k
        return $ GaussLog eg grad


parseEigenData :: GenParser Char st [EigenBLock]
parseEigenData = do
    let eigenLabel = "EIGENVALUES AND  EIGENVECTORS OF CI MATRIX"
        configurations = "NO OF BASIS FUNCTIONS ="
    manyTill anyChar (try $  string configurations)
    cs <- fromInteger `fmap`(spaces >> parseIntegral)
    manyTill anyChar (try $  string eigenLabel)
    manyTill (many1 newline >> parseEigenBLock cs) (try $ string " Final one electron symbolic density matrix:")    
    
parseEigenBLock :: Int -> GenParser Char st EigenBLock
parseEigenBLock n = do
  val <- parseEigenValue
  vec <- parseEigenVector n 
  return $ EigenBLock val vec
      
parseEigenValue :: GenParser Char st Double
parseEigenValue = do
  spaces 
  parenthesis $ spaces >> parseIntegral
  spaces >> string "EIGENVALUE"
  spaces >> parseFloat
  
  
parseEigenVector :: Int -> GenParser Char st [Double]
parseEigenVector n = do 
  pairs <- count n parseConfigVal
  takeLine
  return $ sortPairs pairs  
  
parseConfigVal ::  GenParser Char st (Integer,Double)
parseConfigVal = do
  try spaces <|> (newline >> spaces)
  n <- parenthesis $ spaces >> parseIntegral
  r <- try parseFloat <|> (space >> parseFloat)
  return (n,r)
         
  
-- ======> Parse Gradients <=============

parseLogGrad :: Int -> GenParser Char st  [[Double]]
parseLogGrad numat = do
  let otherState = "Gradient of iOther State \n"
      currentState = "Gradient of iVec State. \n"
  manyTill anyChar (try $ string otherState) 
  otherGrad <- concat `fmap` count numat parseLineNumber 
  manyTill anyChar (try $ string currentState) 
  currentGrad <- concat `fmap` count numat parseLineNumber 
  return [otherGrad,currentGrad]

parseLineNumber :: GenParser Char st [Double]  
parseLineNumber = do
                 xs <- count 3 (spaces >> parseFloat)
                 newline
                 return xs
  
sortPairs :: [(Integer,Double)] -> [Double]
sortPairs = fmap snd . sortBy (compare `on` fst)
  

    
-- =====================> Utilities <===========================    
    

parenthesis p = between (char '(') (char ')')   p

    
-- Parse cardinality field as "N=    Integer"
cardinality :: GenParser Char st Int
cardinality = do
    string "N="
    spaces
    n <- parseIntegral 
    return $ fromIntegral n -- we need just an Int

-- Parse an multi block consisting of n integer values up to six per line
ivalues :: Int -> GenParser Char st [Integer]
ivalues n = 
    if (n `rem` 6 ) == 0 
       then do
           result <- count (n `div` 6) ivaluesline    -- we need just an integral number of lines
           return $ concat result
       else do
           result1 <- count (n `div` 6) ivaluesline   -- we need an integral number of lines 
           result2 <- ivaluesline                     -- plus one more to get the rest of values
           return $ concat [concat result1, result2 ] -- we need to flatten the lists

-- Parse a line of integer values up to six
ivaluesline :: GenParser Char st [Integer]
ivaluesline = do
    idata <- manyTill ivalue newline
    return idata

-- A ivalue consists of one integer after one or more blanks
ivalue :: GenParser Char st Integer
ivalue = do
    simpleSpace
    result <- parseIntegral
    return result

-- A rvalue datablock consists of n rvalues up to six per line
rvalues :: Int -> GenParser Char st [Double]
rvalues n = 
    if (n `rem` 5 ) == 0 
       then do
           result <- count (n `div` 5) rvaluesline
           return $ concat result
       else do
           result1 <- count (n `div` 5) rvaluesline
           result2 <- rvaluesline
           return $ concat [concat result1, result2 ]

-- Parse a line of rvalues up to six
rvaluesline :: GenParser Char st [Double]
rvaluesline = do
    rdata <- manyTill rvalue newline
    return rdata

-- A rvalue consists of one integer after one or more blanks
rvalue :: GenParser Char st Double
rvalue = do
    simpleSpace
    result <- parseFloat
    return result

-- Parse oner or more blanks discarding them
simpleSpace :: GenParser Char st ()
simpleSpace = skipMany (satisfy isSpace)


skipLines :: Int -> GenParser Char st ()    
skipLines n = count n takeLine >> return ()

-- | take whole line
takeLine :: GenParser Char st String
takeLine = manyTill anyChar newline