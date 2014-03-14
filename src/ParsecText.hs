module ParsecText where

import Data.Char
import Text.Parsec
import Text.Parsec.ByteString

-- ========> Internal Imports
import CommonTypes

oneOfStrings :: [String] -> MyParser st String
oneOfStrings listOfStrings = choice $ map (try . string) listOfStrings

-- Parse a text line as a whole into a string
anyLine :: MyParser st String
anyLine = manyTill anyChar newline     -- whatever chars we find till we hit a newline

-- Parser a line containing some string
stringLine :: String -> MyParser st String
stringLine str = do
    spaces
    result <- string str
    manyTill anyChar newline
    return result

blankLines :: MyParser st [String]
blankLines = many1 blankLine

-- Parse a blank line as a whole into a string
blankLine :: MyParser st String
blankLine = manyTill space newline     -- whatever spaces we find till we hit a newline

-- take While symbols are not found without consuming the symbols
takeUntilSymbols :: [Char] -> MyParser st String
takeUntilSymbols xs = many1 $ noneOf xs

manyCharsTillNewLine :: [Char] -> MyParser st String
manyCharsTillNewLine xs = manyTill (oneOf xs) newline

comments, commentsLine :: MyParser st String 
comments     = many1 $ oneOf [' ','*','-','#']
commentsLine = manyTill (oneOf [' ','*','-','#']) newline  
  