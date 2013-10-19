module ParsecText where

import Data.Char
import Text.ParserCombinators.Parsec


oneOfStrings :: [String] -> GenParser Char st String
oneOfStrings listOfStrings = choice $ map (try . string) listOfStrings

-- Parse a text line as a whole into a string
anyLine :: GenParser Char st String
anyLine = manyTill anyChar newline     -- whatever chars we find till we hit a newline

-- Parser a line containing some string
stringLine :: String -> GenParser Char st String
stringLine str = do
    spaces
    result <- string str
    manyTill anyChar newline
    return result

blankLines :: GenParser Char st [String]
blankLines = many1 blankLine

-- Parse a blank line as a whole into a string
blankLine :: GenParser Char st String
blankLine = manyTill space newline     -- whatever spaces we find till we hit a newline
