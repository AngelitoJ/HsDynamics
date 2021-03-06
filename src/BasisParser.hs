
-- Basis Parser, a parsec based parser for quantum chemistry basis files
-- @2013 Angel Alvarez, Felipe Zapata
-- 
--  2013/05/05 We kindly appreciate our mothers efforts in bringing us so far...

-- Limited support for reading STO like basis in nwchem format... 
-- more complex setups will be deployed as needed (do we need uranium right now? (as ))

module BasisParser
    (
        processBasisFile
    ) where

    
import Data.Char
import Data.List as DL

import Text.ParserCombinators.Parsec
import ParsecNumbers
import ParsecText

-- =================================================== Types ===================================================================

type ElementParser st = GenParser Char st Element
type StringParser st  = GenParser Char st String
type AtomLabel = String


data GaussShape = 
      S  Double Double               -- | A S  type Gaussian primitive
    | SP Double Double Double Double -- | A SP type Gaussian primitives contains shared data among S and P types, gaussian functions 
    | P  Double Double               -- | A P  type Gaussian primitive
    | D  Double Double               -- | A D  type Gaussian primitive
    deriving (Show)
    
-- Basic elements living along the lines of a basis file
data Element =
      Header      String              -- | A basis file header spanning just one line
    | Atom        AtomLabel [String] [String] [GaussShape] 
    deriving (Show)

-- =================================================== Main processing functions =========================================================


processBasisFile :: FilePath -> IO ()
processBasisFile fname = do
        putStrLn $ "Processing file:" ++ (show fname) ++ "\n"
        input  <- readFile fname
        case (runParser basisFileParser () fname input) of
             Left err  -> print err
             Right xs  -> mapM_ printElement xs
    
basisFileParser :: GenParser Char st [Element]
basisFileParser = do
    header     <- basisHeader
    references <- many1 basisReferenceSection
    elements   <- parseBasisSection
    eof
    return $ (header:elements) 

-- #  STO-3G  EMSL  Basis Set Exchange Library  5/1/13 6:07 AM
basisHeader :: ElementParser st
basisHeader = do
    char '#'
    many1 space
    basisType <- many1 (noneOf " \n")
    many1 space
    string "EMSL"
    many1 space
    string "Basis Set Exchange Library"
    manyTill anyChar newline
    return $ Header basisType

-- | parse a reference Section as a bunch of references and blanklines
-- | skip those references and blanklines
basisReferenceSection :: GenParser Char st ()
basisReferenceSection = parseReference <|> (skipMany1 blankLine)

-- | parse a Reference block as whole line after a '#'
parseReference :: GenParser Char st ()
parseReference = do
    char '#'
    manyTill anyChar newline
    return ()

-- | A Basis section contains a whole bunch of element and their shape functions
parseBasisSection :: GenParser Char st [Element]
parseBasisSection = do
    string "BASIS \"ao basis\" PRINT"
    newline
    elements <- many1 parseBasisElement
    string "END"
    optional newline
    return $ elements

-- | A element consist of a contraction header and many element's' shape functions sections
parseBasisElement :: ElementParser st
parseBasisElement = do
    (src,dst)      <- atomBasisParameters
    atomShapes     <- many1 atomSection
    let (labels,shapes) = unzip atomShapes
        label           = head labels
        shapeList       = concat shapes
    return $ Atom label src dst shapeList

atomSection :: GenParser Char st (String,[GaussShape])
atomSection = do
    (label,shape) <- atomSectionHeader
    params        <- many1 (shapeParameters shape)
    return $ (label,params)

-- "#BASIS SET: (17s,11p,1d) -> [5s,4p,1d]"
atomBasisParameters :: GenParser Char st ([String],[String])
atomBasisParameters = do
    string "#BASIS SET: "
    sourceG      <- between (char '(') (char ')') ( orbitalIdent `sepBy1` (char ',') )
    string " -> "
    destinationG <- between (char '[') (char ']') ( orbitalIdent `sepBy1` (char ',') )
    newline
    return $ (sourceG, sourceG)

-- | Parse shape parameters of diferent types (with separated or shared exponents)
-- | We have to take care of new shape types if needed...
shapeParameters :: String -> GenParser Char st GaussShape
shapeParameters "S" = do
    (cS,eS) <- parseSingleShape
    return $ S cS eS
shapeParameters "SP" = do                         -- SP types share one parameter
    (cS,cP,eSP) <- parseSharedShape
    return $ SP cS eSP cP eSP
shapeParameters "D" = do
    (cS,eS) <- parseSingleShape
    return $ D cS eS
shapeParameters sh = do
    unexpected $ "Gaussian Shape type " ++ sh ++ " not currently supported, (just tell us)"
    return $ S 0.0 0.0

-- | Parse two separate parameters (single types S,P,D etc...)
-- "  25180.1000000              0.0018330        "
parseSingleShape :: GenParser Char st (Double,Double)
parseSingleShape = do
    many1 space
    eA            <- realNumber
    many1 space
    cA            <- realNumber
    manyTill space newline
    return (cA,eA)

-- | Parse three parameters (for shared types like SP)
-- "      3.1864900             -0.2518300             -0.0142990        "
parseSharedShape :: GenParser Char st (Double,Double,Double)
parseSharedShape = do
    many1 space
    eAB           <- realNumber
    many1 space
    cA            <- realNumber
    many1 space
    cB            <- realNumber
    manyTill space newline
    return (cA,cB,eAB)

-- "Cl    SP"
atomSectionHeader :: GenParser Char st (String,String)
atomSectionHeader = do
    label <- oneOfStrings elementStrings  -- Try every "known" element in order
    skipMany1 space
    shape <- oneOfStrings $ lowerAndUpper orbitalStrings  -- Try every "known" shape in order
    newline
    return $ (label,shape)
    

orbitalIdent :: StringParser st 
orbitalIdent = do
    number <- intNumber
    shape  <- oneOfStrings orbitalStrings
    return $ (show number) ++ shape

-- =================================================== Utility functions ===================================================================

-- | convert a list o lowercase estring into a list with both lower and uppercase elements
lowerAndUpper :: [String] -> [String]
lowerAndUpper = (foldl (\a x -> x:(map toUpper x):a ) []).(reverse)

-- | Known list of molecular orbital types. People sugest there are more...
orbitalStrings :: [String]
orbitalStrings =  ["d","f,","p","sp","s"]

-- | Western parties publicy available elements list. Commies seemed to know more elements but all knowledge got lost after the USSR passed away
elementStrings :: [String]
elementStrings = ["Be","Ca","Cl","C","H","K","Li","Na","N","Mg","O"]

-- | Convert beautifull pure functional data into nasty efects...
-- | Keep elegance whenever posible
printElement :: Element -> IO ()
printElement (Header   header) = putStrLn $ "Basis Set: " ++ header
printElement (Atom     label src dst shapes) = do
    putStrLn $ "Atom "++ label ++ " " ++ (show src) ++ " -> " ++ (show dst)
    mapM_ printShape shapes
    where
        printShape sh = putStrLn $ "\t" ++ ( show sh)
