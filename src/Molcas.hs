-- ------------------------------------------------------------------------------------
-- Molcas 2010 output file parser
-- @2013 The ResMol Group
-- Authors: Angel Alvarez, Felipe Zapata
--
--  2013-01 Initial design

-- Molcas:
-- Another bunch of ancient and (somewhat) heroic Fortran routines dating back from the 
-- mid-eighties to nowadays, yet a bit more structured that Gaus** routines that look 
-- even older... Both share the batch model of execution where multiple programs run in 
-- sequence to acomplish some calculation and (to some extent) they exercise common routines.
-- 
-- The output file:
-- The format seems a bit chaotic (even for chemistry?) and poorly suited for machine parsing. 
-- Yet lot of people try to make a life upon this messy outputs so it deserve a good parser.
-- 
-- YMMV comes to be the norm when expecting where data will be and why, actually some info 
-- moves back and forth between runs even to unexpected places, upon many different facts
-- not always related to alternative code branchs. 

-- Such a mess makes you feel shocked that this programming culture still achieves so many 
-- things. Fortran is no longer the king in the compiler arena, neither is the language of choice.
-- Today a plethora of developer tecnologies achieve impressive performance yet leveraging good 
-- software engineering and parallel hardware is a commodity. Surely there are better and more 
-- efficient ways nowadays to building number crunching and theoretic software so here we are.
-- 
-- Even as is, this is a piece of impressive knowledge and hard work, we honor that. 
-- So please take no ofense. Hopefully this work will somehow allow people easy data parsing of 
-- this gigant in chemistry.
--
--                     ....We promise, no chemies were harmed during this development.
-- ------------------------------------------------------------------------------------

module Molcas where

import Control.Applicative ((<$>),(<*>),(*>),(<*))
import Control.Exception (assert)
import Data.Char
import qualified Data.ByteString.Char8 as C
import qualified Data.List.Split as DLS
import Text.Parsec
import Text.Parsec.ByteString

-- ======> Internal imports <=========
import CommonTypes
import ParsecNumbers
import ParsecText

-- ==========================================================

-- We try to parse and convert unformatted text data into more manageable chunks with gracefull degradation
-- 
-- 1) Highly formatted (i.e. Molecular gradients) structures upon tuples and primtive types.
-- 2) Medium formatted data (i. e. Alaska data) where we store data we havent work it out yet.
-- 3) Unformatted data, we contruct text chunks out of text lines, yet we try to record where we found it 
-- 
-- As soon as a project of Us needs new data we make proper parsing combinators and types to support it.


-- Public Interface ----------------------------------------------------------------------------------


-- Create a defualt Molcas Parser state 
defaultParserState :: MolState 
defaultParserState = MolState {
                              blockLevel   = []  
                            , sectionLevel = []    
                        }

-- Parse a given file path containing a ouput Molcas file and display it
processMolcasOutputFile :: FilePath -> IO ()
processMolcasOutputFile fname = do
        putStrLn $ "Processing file:" ++ (show fname) ++ "\n"
        input  <- C.readFile fname
        case (runParser molcasOutParser defaultParserState fname input) of
             Left err  -> print err
             Right xs  -> mapM_ print xs

-- Parse a given file path pointing to a Molcas output file 
parseMolcasOutputFile :: FilePath -> IO (Either ParseError [MolBlock])
parseMolcasOutputFile fname = do
        input  <- C.readFile fname
        return $ runParser molcasOutParser defaultParserState fname input
        
pruebaParser :: Show a => FilePath -> MyParser MolState a -> IO ()
pruebaParser  fileName parser = do
  input <- C.readFile fileName
  case runParser parser defaultParserState fileName input of
       Left msg -> print msg
       Right xs -> print xs

-- ------------------------------------------------------------------------------------------------------------
--  Top level parsers, one to input files (just a hack), one to Molcas output files 
-- ------------------------------------------------------------------------------------------------------------

-- Quick Hack to parse the Molcas Input file format and gather a parameter
molcasInpParser :: MyParser MolState (String,Int,String)
molcasInpParser = do
    firstPart <- manyTill anyChar (try (string "SALA"))
    newline
    spaces
    root <- intNumber
    spaces
    lastPart <- manyTill anyChar (try eof)
    return $ (firstPart,root,lastPart)

-- Parsec parser for Molcas Output format
molcasOutParser :: MyParser MolState [MolBlock]
molcasOutParser = do
    modifyState (\x -> x { blockLevel = parsers } ) -- Set the Block parsers according to main part of a Molcas file
    elements <- many1 molcasBlock                   -- Get at least one of the known Molcas Blocks
    return $ elements
    where
        parsers = [ 
                    molcasBlank                       -- Who needs blanks? we don't
                    , molcasLicense                   -- Is this a Molcas license block?
                    , molcasLogo                      -- or the Molcas "Logo"?
                    , molcasCopyRight                 -- or who rules the world?
                    , molcasProject                   -- or project data?
                    , molcasAuto                      -- the auto module
                    , molcasWhatever                  -- You guess what? Actually I just dont care
                    ]

 -- ---------------------------------------------------------------
 -- Molcas Input Parsers
 -- 
-- -----------------------------------------------------------------
parserInputMolcasQM :: FilePath -> MyParser MolState [AtomQM] -> IO [AtomQM]
parserInputMolcasQM  fileName parser = do
  input <- C.readFile fileName
  case runParser parser defaultParserState fileName input of
       Left msg -> error $ show msg
       Right xs -> return xs

-- | Molcas .input File Parser       
parseMolcasInputFile :: FilePath -> IO [MolcasInput String]
parseMolcasInputFile fname = do
  input  <- C.readFile fname
  case runParser parserMolcasInput defaultParserState fname input of
       Left msg -> error $ show msg
       Right xs -> return xs       

parserMolcasInput :: MyParser MolState [MolcasInput String]
parserMolcasInput = many1 parseInputSection 
  

parseInputSection :: MyParser MolState (MolcasInput String)
parseInputSection =    try parserInputCommand
                   <|> try parserGateway
                   <|> try parserSeward
                   <|> try parserCasscf
                   <|> try parserMclr
                   <|> try parserAlaska
                   <?> "I don't know which Molcas Module are you talking about, is it Swedish?"  

parserInputCommand :: MyParser MolState (MolcasInput String )
parserInputCommand = do
                     string ">>"
                     Command <$> anyLine

parserGateway :: MyParser MolState (MolcasInput String )
parserGateway = do
                manyTill anyChar $ char '&'
                oneOf ['G','g'] >> anyLine
                Gateway <$> (takeUntilSymbols ['&'])

parserSeward :: MyParser MolState (MolcasInput String )
parserSeward = do
               manyTill anyChar $ char '&'
               oneOf ['S','s'] >> anyLine
               Seward <$> (takeUntilSymbols ['&'])

parserCasscf :: MyParser MolState (MolcasInput String )
parserCasscf = do
               manyTill anyChar $ char '&'
               oneOf ['R','r'] >> anyLine
               ini  <- manyTill anyChar (try (string "rlxroot") <|> try (string "Rlxroot") <?> "expecting Relax Root keyword" )
               spaces >> char '=' >> spaces
               n    <- intNumber
               rest <- takeUntilSymbols ['&']
               return $ RasSCF n ini rest

parserMclr :: MyParser MolState (MolcasInput String )
parserMclr = do
             manyTill anyChar $ char '&'
             oneOf ['M','m'] >> anyLine
             MCLR <$> (takeUntilSymbols ['&'])
             

parserAlaska :: MyParser MolState (MolcasInput String )
parserAlaska = do
               manyTill anyChar $ char '&'
               (oneOf ['A','a']) >> anyLine
               return $ Alaska ""
--                manyTill anyChar $ try (oneOf ['&','>']) <|> try eof  


isInputRasscf :: (MolcasInput String) -> Bool
isInputRasscf x = case x of
                       RasSCF _ _ _ -> True
                       otherwise    -> False

parserGatewayQM :: Int -> MyParser MolState [AtomQM]
parserGatewayQM n = try $ parserAtoms n

parserAtoms:: Int -> MyParser MolState [AtomQM]
parserAtoms numat = do 
    manyTill anyChar (  try (string "&Gate")
                    <|> try (string "&gate")
                    <|> try (string "GATE" )
                    <?> "I was looking for the Gateway keyword" )
    anyLine                       
    count numat parserAtom 
    
parserAtom :: MyParser MolState AtomQM
parserAtom =        try parserAtomMM
                <|> try parserAtomQM
                
parserAtomMM :: MyParser MolState AtomQM                 
parserAtomMM = do    
    spaces >> (try (string "Basis set\n") <|> try (string "basis set\n"))   -- header 
    spaces >> anyChar >> string "...... / MM\n"                             -- MM linker 
    label <- spaces >> manyTill anyChar digit                               -- atomic label
    manyTill digit space                                                    -- rest of the label numbering
    xyz   <- count 3 (spaces >> realNumber)                                 -- Cartesian coordinates
    spaces >> try (string "Angstrom" <|> string "angstrom")
    spaces >> try (string "Charge" <|> string "charge") >> spaces >> char '=' >> spaces >> realNumber
    spaces >> try (string "End of Basis\n")
    return $ AtomQM label xyz
    
    
parserAtomQM :: MyParser MolState AtomQM                 
parserAtomQM = do
     spaces >> try (string "Basis set\n" <|> string "basis set\n")
     anyLine  -- basis Line
     label <- spaces >> manyTill anyChar digit
     manyTill digit space
     xyz   <- count 3 (spaces >> realNumber)
     manyTill anyChar newline
     spaces >> try (string "End of Basis\n")
     return $ AtomQM label xyz


     
 -- ------------------------------------------------------------------------------------------------------------
--  Molcas Output Block level parsers:
--  one to every know main block of a output Molcas file, ie: License,Project info and so on
-- ------------------------------------------------------------------------------------------------------------


-- Parse a Molcas block from the (user state provided) list of known block level parsers 
molcasBlock :: BlockParser
molcasBlock = do
    userState <- getState                 -- get block level parsers from the user state
    choice $ blockLevel userState           -- try the top level parsers
 
-- Try to parse as many blank text lines as it can, and collapse them into a single TextData MolBlock 
molcasBlank :: BlockParser
molcasBlank = try $ do
    many1 newline                   -- take easy man, we dont judge what we find...
    return $ TextData "[BLANK]"

-- Parse while text lines into MolBlocks 
molcasWhatever :: BlockParser
molcasWhatever = do
    iFound <- anyLine                   -- take easy man, we dont judge what we find...
    return $ TextData iFound

-- Parse Molcas License Info
molcasLicense :: BlockParser
molcasLicense = try $ do
    string "   License is going to expire in "
    days      <- intNumber
    string " days  on "
    date      <- manyTill anyChar newline
    string "   This copy of MOLCAS is licensed to "
    user      <- manyTill anyChar newline
    let licenseData = "Expires in " ++ date ++ " (" ++ (show days) ++ "), Licensed to: " ++ user
    return $ LicenseData licenseData

-- Parse Molcas Logo 
molcasLogo :: BlockParser
molcasLogo = try $ do
    spaces
    string "^^^^^"
    spaces
    string "M O L C A S"
    spaces
    spaces
    string "^^^^^^^"
    spaces
    string "version"
    space
    version       <- versionString -- 7.6
    string " patchlevel "
    patch         <- manyTill digit space
    spaces
    count 35 anyLine              -- its a big logo!
    let logoData = "Version " ++ version ++ " patchlevel " ++ patch
    return $ LogoData logoData    -- return as single string

-- Try parse a Molcas copyright text block, return only university info and year
molcasCopyRight:: BlockParser
molcasCopyRight = try $ do
    stringLine "Copyright, all rights, reserved:"
    stringLine "Permission is hereby granted to use"
    stringLine "but not to reproduce or distribute any part of this"
    stringLine "program. The use is restricted to research purposes only."
    message <- stringLine "Lund University Sweden, 2010."
    stringLine "For the author list and the recommended citation consult section 1.5 of the MOLCAS user's guide."
    return $ CopyRightData message

    
-- Try parse the input file data from the auto module ouput
molcasInput :: BlockParser
molcasInput = try $ do
    string "++ ---------   Input file   ---------"
    contents <- manyTill anyLine (try (string "-- ----------------------------------"))
    newline
    return $ InputData (concat contents)

--  Try to parse Project data returning a proper ProjectData block
molcasProject:: BlockParser
molcasProject = try $ do
    string "   -------------------------------------------------------------------"
    newline
    string "  |                                                                   "
    newline
    string "  |   Project         = "
    name                   <- manyTill anyChar (try newline) -- "gradR4"
    string "  |   Submitted from  = "
    submitted              <- manyTill anyChar (try newline) -- "/scratch/cris/gradR4.3173"
    string "  |   Scratch area    = "
    scratch                <- manyTill anyChar (try newline) -- "/scratch/cris/gradR4.3173/"
    string "  |   Save outputs to = "
    output                 <- manyTill anyChar (try newline) -- "/home/cris/MOL_MOT/luisma_iminium_out/conform/conf_S0_dul/gradientes/R4"
    string "  |                                                                   "
    newline
    string "  |   Scratch area is NOT empty"
    newline
    string "  |                                                                   "
    newline
    string "  |                                                                   "
    newline
    string "   -------------------------------------------------------------------"
    newline
    return $ ProjectData name submitted scratch output

-- Try to parse the auto module output, store runtime statistics and include output inside module runs
molcasAuto :: BlockParser
molcasAuto = try $ do
    (name, whenStart)             <- moduleStart                           -- The module started
    blankLine
    input                         <- molcasInput                           -- A copy of input file is provided
    blankLine
    workdir                       <- option "" workDir                     -- One or more workdir sentences can appear here (and there)
    modules                       <- many1 molcasModule                    -- A bunch of modules, come ahead...
    run                           <- autoStatus                            -- Dont look at me!, im just parsing your results...
    (name1, whenStop,returnCode)  <- moduleStop                            -- The module stopped, hopefully
    duration                      <- option "" $ moduleSpent name
    return $ ModuleAuto name (whenStart ++ " " ++ whenStop ++ " " ++ duration ++ " RC: " ++ returnCode) modules


-- Parse a normal (non auto) Molcas Module, and check for suitable sections upon checking module name
-- This parsing is crazy, Molcas modules seem to have duplicate branchs of code printing out the same
-- info yet not really the same...
molcasModule :: MyParser MolState ModuleData
molcasModule = try $ do
    (name, whenStart)                 <- moduleStart                    -- The module started
    items                             <- manyTill (sectionParser name) $ lookAhead $ try moduleStop 
    (moduleName, endTime, returnCode) <- moduleStop
    duration                          <- option "" $ moduleSpent name   -- sometimes it just pops up, sometimes not 
    workdirReset                      <- option "" workDir              -- Who is messing with the work dir, uhh?
    autospent                         <- option "" $ moduleSpent "auto" -- "Shit happens" - Forrest Gump
    return $ ModuleData name ("Started: " ++ whenStart ++ " Stopped (" ++ returnCode ++ ") " ++ endTime ++ " " ++ duration ) items

-- "--- Start Module: seward at Sun Dec  2 22:41:44 2012"
moduleStart :: MyParser MolState (String,String)
moduleStart = do
    optional ( try (string "*** " >> newline))                       -- Depends upon weather on Sweden, the module and the daily fellow 
    string "--- Start Module:"                                       -- Here we go, fasten your seat belts
    skipMany space
    moduleName <- many1 alphaNum
    skipMany space 
    string "at"
    content <- anyLine
    return $ (moduleName, content )

-- "--- Stop Module:  gateway at Sun Dec  2 22:41:43 2012 /rc=0 ---" 
moduleStop :: MyParser MolState (String, String, String)
moduleStop = do
    string "--- Stop Module:"
    skipMany space
    moduleName <- many1 alphaNum
    skipMany space 
    string "at"
    endTime  <- manyTill anyChar (string "/rc=")
    returnCode <- manyTill anyChar (string "---" >> newline)
    return $ (moduleName, endTime, returnCode)

-- Try to guess how this run of Molcas finished, return a (somewhat) descriptive message
autoStatus :: MyParser MolState String
autoStatus = do
    happyLanding 
    <|> convergenceProblem                                            -- You fool, that's pretty big for us
    <|> unknownProblem                                                -- Cosmic Rays
    <?> "Really hard to guess this!"
    where
        happyLanding :: MyParser MolState String
        happyLanding = try $ do
            blankLine
            string "     Happy landing!"                               -- Normal case you are lucky such a beast didnt kill your cat 
            newline
            blankLine
            return "Happy Landing"
        convergenceProblem :: MyParser MolState String
        convergenceProblem = try $ do
            string "          "
            newline
            string "    Convergence problem!"                          -- Error Convergence problem, whatever that means 
            newline
            string "          "
            newline
            count 8 anyLine                                            -- Molcas Gurus wont tell you, man you know you did wrong..
            string "Non-zero return code - check program input/output" -- Clever clue
            newline
            return "Convergence Problem"
        unknownProblem :: MyParser MolState String
        unknownProblem = try $ do
            many1 anyLine                                              -- Who knows, just ask Roland...
            return "Unknown Problem"

-- Parse WorkDir rset line
-- This message  is know to appear in many unrelated parts of the output file...
workDir :: MyParser MolState String
workDir = do
    lines <- many1 (string "WorkDir " >> anyLine)
    return $ head lines

-- Parse a simple version string as [0-9]'.'[0-9]
versionString :: MyParser MolState String
versionString = do
    major <- digit
    char '.'
    minor <- digit
    let version = major : '.' : [minor]
    return $ version
    

moduleSpent :: String -> MyParser MolState String
moduleSpent name = try $ do
    string "--- Module "
    string name
    string " spent "
    spent <- anyLine
    return spent

-- ------------------------------------------------------------------------------------------------------------
--  Molcas Output Section level parsers:
--  one to every know sections from modules data, i.e. Alaska Molecular Gradients or Rasscf Energies
-- ------------------------------------------------------------------------------------------------------------

sectionParser :: String -> SectionParser
sectionParser name = case name of
    "alaska" -> do
                 sectionMolecularGradients                  -- Molecular Gradient printout 
                <|> sectionsWhatever name                  -- rest of Alaska sections
                <?> "alaska me cachis!"                    -- You know what I mean, man
    "rasscf" -> do
                sectionCiCoefficients                      -- rasSCF Wave function printout (occupation of active orbitals, and spin coupling of open shells )
                <|> sectionRASSCFEnergies                  -- Final rasSCF Energies
                <|> sectionsWhatever name                  -- rest of RASSCF sections
                <?> "rasscf me cachis!"                    -- Should we phone Roland one more time?
    other    -> do
                sectionsWhatever other                     -- We dont have a section parser for this module sections...

-- Parse sections lines into a suitable generic ModuleData 
sectionsWhatever :: String -> SectionParser
sectionsWhatever "alaska" = do
    iFound <- anyLine
    return $ AlaskaData iFound
sectionsWhatever "rasscf" = do
    iFound <- anyLine
    return $ RasSCFData iFound
sectionsWhatever name = do
    iFound <- anyLine
    return $ WhateverData $ "[" ++ name ++ "] " ++ iFound    -- Well, I dont know what is it


-- Try to parse Molecular gradients from alaska Molcas modules
sectionMolecularGradients :: SectionParser
sectionMolecularGradients = try $ do
    string " **************************************************"
    newline
    string " *                                                *"
    newline
    string " *              Molecular gradients               *"
    newline
    string " *                                                *"
    newline
    string " **************************************************"
    newline
    blankLine
    string "                Irreducible representation:"
    spaces
    representation <- oneOfStrings ["a","s1"]
    spaces
    elements       <- many1 molecularGradient
    let n = (`div`3) $ length elements
    espf           <- option [] $ try $ parserGradESPF n
    return $ AlaskaMolecularGradients representation elements espf


-- Try parse molecular gradients items from alaska
-- "                C1       x                 0.6259617E-01"
molecularGradient:: MyParser MolState (String,Char,Double)
molecularGradient = try $ do
    spaces
    atom           <- many1 alphaNum
    spaces
    coord          <- oneOf "xyz"
    spaces
    value          <- realNumber
    newline
    return $ (atom,coord,value)
    
parserGradESPF :: Int -> MyParser MolState [[Double]]
parserGradESPF n = do
    manyTill anyChar (try $ string "After ESPF")  
    count 2 (anyLine)
    count n parserLineXYZ
   
parserLineXYZ :: MyParser MolState [Double]
parserLineXYZ = do 
             spaces >> intNumber 
             (count 3 $ spaces >> realNumber)

-- Try parse Ci root data from rasscf Molcas module
sectionCiCoefficients :: SectionParser
sectionCiCoefficients = try $ do
    skipMany1 space
    string "printout of CI-coefficients larger than"
    skipMany1 space
    threshold <- realNumber                   --"0.05" 
    skipMany1 space
    string "for root"
    skipMany1 space
    root <- intNumber
    newline
    skipMany1 space
    string "energy="
    skipMany1 space
    energy <- realNumber
    newline
    skipMany1 space
    string "conf/sym"
    skipMany1 space
    mask <- manyTill (char '1') space -- 11111111
    skipMany1 space
    string "Coeff"
    skipMany1 space
    string "Weight"
    newline
    elements <- many1 coefficientElement
    return $ RasSCFCI root energy elements

sectionRASSCFEnergies :: SectionParser
sectionRASSCFEnergies = try $ do
    string "      Final state energy(ies):"
    newline
    string "      ------------------------"
    newline
    anyLine
    elements <- many1 rootEnergy
    return $ RasSCFRE elements

-- Parse rasscf root Energy
-- ::    RASSCF root number  2 Total energy =       -821.43755900                                                          
rootEnergy :: MyParser MolState (Int,Double)
rootEnergy = try $ do
    string "::    RASSCF root number"
    spaces
    root <- intNumber
    string " Total energy ="
    spaces
    energy <- realNumber
    anyLine
    return $ (root,energy)

-- Parse rasscf CI root info "    1  22220000  -0.91395 0.83531"
coefficientElement :: MyParser MolState (Int, String, Double, Double)
coefficientElement = try $ do
    skipMany1 space
    confsym <- intNumber
    skipMany1 space
    flags <- many1 alphaNum
    skipMany1 space
    coeff <- realNumber
    skipMany1 space
    weight <- realNumber
    newline
    return (confsym, flags, coeff, weight)

-- -------------------------------------------------------------------------------
-- 
-- Filtering combinators ans Show Instances 
-- 
-- -------------------------------------------------------------------------------

-- is this Block, the Molcas auto module?
isAuto :: MolBlock -> Bool
isAuto (ModuleAuto _n _s _modules) = True
isAuto _ = False

-- Test if given Module have a specified name
isModule :: String -> ModuleData -> Bool
isModule name (ModuleData n _s _sections) = n == name

--  Combine two module tests into one OR'ing the results of both test
orModule :: (ModuleData -> Bool) -> (ModuleData -> Bool) -> (ModuleData -> Bool)
orModule a b = \x -> (a x) || (b x)

infixr 2 `orModule`

-- Test for desired Module
isAlaska, isRASSCF :: ModuleData -> Bool
isAlaska = isModule "alaska"
isRASSCF = isModule "rasscf"

-- test for desired SectionData
isCICoeffs, isREnergies, isMolGrads :: SectionData -> Bool
isCICoeffs (RasSCFCI _r _e _ci) = True
isCICoeffs _ = False
isREnergies (RasSCFRE _roots) = True
isREnergies _ = False
isMolGrads (AlaskaMolecularGradients _r _grads _gradESPF) = True 
isMolGrads _ = False


