



module GenericParser 
    (
        processGenericFile
    )where


import Text.ParserCombinators.Parsec
import ParsecNumbers


data Element =
      RData Double
    | IData Int
    | SData String
    | CData Char
    | Space Int
    
data FileLine = FileLine Int String [Element]

type StringParser st  = GenParser Char st String
type LineParser st    = GenParser Char st FileLine
type ElementParser st = GenParser Char st Element

-- instance Show FileLine where
--     show = showLine
-- 
-- instance Show Element where
--     show = showElement
-- 
-- showLine:: FileLine -> String
-- showLine (FileLine l d elements) = "[" ++ (show l) ++ "] {" ++ d ++ "} ->" ++ (concat $ map show elements) ++ "<-"
-- 
-- showElement:: Element -> String
-- showElement (RData _x) = "R"
-- showElement (IData _x) = "I"
-- showElement (SData _x) = "S"
-- showElement (CData _x) = "C"
-- showElement (Space 0)   = "*"
-- showElement (Space 1)   = "_"
-- showElement (Space x)   = "_" ++ show(x) ++ "_"

processGenericFile :: FilePath -> IO ()
processGenericFile file = do
        putStrLn $ show file
        result <- parseFile file
        case result of
             Left err  -> print err
             Right xs  -> mapM_ dumpLine xs

dumpLine :: FileLine -> IO ()
dumpLine l@(FileLine i d e) = do
    putStrLn $ (show i) ++ " {" ++ d ++ "} " ++ linetrace ++ "     " ++ linecontents
    where
        linetrace = concat $ map traceElement e
        linecontents = concat $ map dumpElement e
        traceElement:: Element -> String
        traceElement (RData _x) = "R"
        traceElement (IData _x) = "I"
        traceElement (SData _x) = "S"
        traceElement (CData _x) = "C"
        traceElement (Space 0)   = "*"
        traceElement (Space 1)   = "_"
        traceElement (Space x)   = "_" ++ show(x) ++ "_"
        dumpElement:: Element -> String
        dumpElement (RData x) = show x
        dumpElement (IData x) = show x
        dumpElement (SData x) = x
        dumpElement (CData x) = [x]
        dumpElement (Space 0) = ""
        dumpElement (Space 1) = " "
        dumpElement (Space x) = take x $ repeat ' '


parseFile :: String -> IO ( Either ParseError [FileLine])
-- parseFile = parseFromFile parseGenericFile 
-- parseFromFile :: Parser a -> String -> IO (Either ParseError a)
parseFile fname = do
    input <- readFile fname
    return $ runParser parseGenericFile () fname input



parseGenericFile :: GenParser Char st [FileLine]
parseGenericFile = do
    elements <- many parseFileLine
    eof
    return $ elements

parseFileLine :: LineParser st
parseFileLine = do
        parseNewLine
    <|> parseBlankLine
    <|> parseLine
    <|> parseAnyLine
    <|> parseAnyRest
    <?> "kaka"

parseNewLine :: LineParser st
parseNewLine = try $ do
    pos <- getPosition
    newline
    return $ FileLine (sourceLine pos) "Newline" [Space 0]

parseBlankLine :: LineParser st
parseBlankLine = try $ do
    pos <- getPosition
    elements <- space `manyTill` newline
    return $ FileLine (sourceLine pos) "BlankLine" [Space (length elements)]

parseLine :: LineParser st
parseLine = try $ do
    pos <- getPosition
    elements <- parseElement `manyTill` newline
    return $ FileLine (sourceLine pos) "Line" elements

parseAnyLine :: LineParser st
parseAnyLine = try $ do
    pos <- getPosition
    elements <- anyChar `manyTill` newline
    return $ FileLine (sourceLine pos) "anyLine" [SData elements]

parseAnyRest :: LineParser st
parseAnyRest = try $ do
    pos <- getPosition
    elements <- many1 anyChar
    return $ FileLine (sourceLine pos) "anyRest" [SData elements]

parseElement :: ElementParser st
parseElement = do
    parseSpace
    <|> parseNum
    <|> parseWord

parseSpace :: ElementParser st
parseSpace = do
    content <- many1 (char ' ')
    return $ Space (length content)

parseNum :: ElementParser st
parseNum = do
        parseReal
    <|> parseNat

parseNat :: ElementParser st
parseNat = try $ do
    number <- intNumber
    return $ IData number

parseReal :: ElementParser st
parseReal = try $ do
    number <- realNumber
    return $ RData number

-- Parse a single char or a string of printable caracters
parseWord :: ElementParser st
parseWord = try $ do
    content <- many1 (noneOf " \n")
    return $ charOrString content
    where
        charOrString x  
            | (length x > 1) = SData x          -- Many matches so we return a string
            | otherwise      = CData (head x)   -- Single match so we return a char
        
-- parseChar :: ElementParser st
-- parseChar = try $ do
--     c <- anyChar
--     return $ CData c
