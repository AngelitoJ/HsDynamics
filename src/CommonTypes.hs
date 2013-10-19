{-# Language DeriveFunctor #-}


module CommonTypes where

import Data.Array.Repa as R hiding ((++))
import Data.Complex
import Data.Functor.Identity
import Data.List as DL
import qualified Data.List.Split as DLS
import qualified Data.Map as M
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector as V
import Text.Parsec


-- ========================> Types for Parsing <======================================

data BlockType = GauBlock | MolBlock deriving Show     

data EigenBLock = EigenBLock { getRootEnergy :: !Double 
                              ,getRootCoeff ::  ![Double] }  deriving Show

data GauBlock = 
      IGauBlock Label Int [Integer]     -- A GauBlock containing one or more integer values
    | RGauBlock Label Int [Double]      -- A GauBlock containing one or more real values
    | TGauBlock String                  -- A GauBlock containing just unformatted text

    
data GaussLog = GaussLog [EigenBLock] [[Double]] deriving Show
    
data InitialDynamics = InitialDynamics {
                                        getInitialState :: Singlet
                                       ,getTime :: Double
                                       ,getdt   :: Double
                                       ,getExtForceMod :: Double 
                                       ,getForceAnchor :: [Int]
                                       ,getTheoryLevel :: String
                                       }    deriving Show

                                       
-- ==========> Molcas Types <=========                                      
type BlockParser   = MyParser MolState MolBlock
type SectionParser = MyParser MolState SectionData 

-- A parsecT based parser carring state of type "a" and returning data of type "b"
type MyParser a b  = ParsecT [Char] a Identity b

-- Pretty simple types right now, but will have type parameters in the near future
-- Molblocks represent main Molcas format blocks including the auto module as a whole
data MolBlock =
      LicenseData String                                     -- LicenseData contains license information 
    | LogoData    String                                     -- LogoData contains version info 
    | CopyRightData String                                   -- Copy... contains law related things
    | InputData  String                                      -- Original input file text used for this run
    | ProjectData String String String String                -- Project details, submission, scratch dir etc...
    | ModuleAuto Name String [ModuleData]                    -- Molcas auto module and all moules it ran
    | TextData    String                                     -- just unformatted text as we got it 

-- Most Molcas modules run from the auto module since the latter control flow of execution
data ModuleData = 
     ModuleData Name String [SectionData]                    -- Molcas modules data stored in sections 

-- Every molcas module produces some data results, some of them have diferenciated types for better exploitation...
data SectionData = 
      WhateverData String                                    -- Data not really parsed and left as unformatted data
    | GatewayData String
    | SewardData String
    | RasSCFData String                                      -- RASSCF unformatted output
    | RasSCFCI Int Double [(Int, String, Double, Double)]    -- RASSCF CI Roots info
    | RasSCFRE [(Int,Double)]                                -- RASSCF Energies
    | AlaskaData String                                      -- Alaska unformatted output
    | AlaskaMolecularGradients String [(String,Char,Double)] -- Alaska Molecular Gradients

-- Common state record to all Molcas parsers
data MolState = MolState 
    { 
          blockLevel   :: [BlockParser]                       -- Block level Molcas blocks
        , sectionLevel :: [SectionParser]                     -- Section Level Molcas Blocks
    }
    
data Job = Gaussian TheoryLevel | Interpolation | Molcas FilePath | Quadratic | HaskellAbInitio deriving Show 


type Anchor = [Int]
type Args = String
type Command = String
type Coordinates = Array U DIM1 Double
type DT = Double -- step of integration
type Energy = Double
type Label = String
type Masses = Array U DIM1 Double
type MatrixCmplx = [Complex Double]
type Name = String
type Project = String
type ReactionCoordinate = M.Map Int (Int,Internals)
type Roots       = Int    -- Number of roots of the system
type TheoryLevel = String
type Temperature = Double
type TotalEnergy = Double
type Time = Double
-- ================> Instances <=================
instance Show MolBlock where
    show (TextData s)             = "Text    : " ++ s
    show (LicenseData s)          = "License : " ++ s
    show (LogoData s)             = "Logo    : " ++ s
    show (CopyRightData s)        = "Copy Right : " ++ s
    show (InputData s)            = "Input   : " ++ s
    show (ProjectData p f s x)    = "Project info:\n\tName: " ++ p ++ "\n\tSubmitted: " ++ f ++ "\n\tScratch: " ++ s ++ "\n\tOuputs: " ++ x
    show (ModuleAuto n s modules) = "Module : " ++ n ++ " " ++ s ++ "\n" ++ concat (fmap show modules)

instance Show ModuleData where
    show (ModuleData n s items)   = "Module : " ++ n ++ " " ++ s ++ "\n" ++ concat (fmap show items)

instance Show SectionData where
    show (WhateverData s)                   = "Module Data: " ++ (show s) ++ "\n"
    show (RasSCFData s)                     = "RasSCF Data: " ++ (show s) ++ "\n"
    show (RasSCFCI r e ci)                  = "RasSCF CI-Coefficients for root: " ++ (show r) ++ " Energy: " ++ (show e) ++ " info:" ++ (show ci) ++ "\n"
    show (RasSCFRE roots)                   = "RasSCF Final Energies:" ++ (show roots)++"\n"
    show (AlaskaData s)                     = "Alaska Data: " ++ (show s) ++ "\n"
    show (AlaskaMolecularGradients r grads) = "Alaska Molecular Gradients [" ++ r ++ "]:" ++ (showMolGrads grads) ++ "\n"

-- Make a Alaska Molecular Gradients listing suitable for humans...
showMolGrads :: [(String,Char,Double)] -> String
showMolGrads list = concatMap showGrad $ chunkedlist
    where
        chunkedlist = DLS.chunksOf 3 list
        showGrad [(a1,x,vx),(a2,y,vy),(a3,z,vz)] = "\n\tAtom " ++ a1
            ++ "  " ++ [x] ++ " = " ++ (show vx)
            ++ ", " ++ [y] ++ " = " ++ (show vy)
            ++ ", " ++ [z] ++ " = " ++ (show vz)


-- ==================> Internal Coordinates <============
data EnergyDerivatives = EnergyDerivatives !Energy !Grad !Hess deriving Show

data FrankCondon = FC !Internals !Connections deriving Show

data  InternalCoord = Bond !Int !Int | Angle !Int !Int !Int | Dihedral !Int !Int ! Int !Int deriving Show                          

type Connections = V.Vector InternalCoord
type Grad = Array U DIM1 Double
type Hess = Array U DIM2 Double
type Internals = VU.Vector Double



-- =================> Types for Dynamics <===============

data Singlet = S0 | S1 | S2 | S3 | S4 deriving (Show,Enum,Eq,Read)

data Triplet = T1 | T2 | T3 | T4 | T5 deriving (Show,Enum,Eq,Read)


-- | Data representing some values associated to a molecule
data Molecule = Molecule {
                getCoord   :: !(R.Array R.U R.DIM1 Double)
               ,getVel     :: !(R.Array R.U R.DIM1 Double)
               ,getForce   :: !(R.Array R.U R.DIM1 Double)
               ,getMass    :: !(R.Array R.U R.DIM1 Double)
               ,getAtoms   :: ![String]
               ,getEnergy  :: ![[Double]]
               ,getDervEn  :: ![EnergyDerivatives]
               ,getElecSt  ::  Either Singlet Triplet
               ,getFCStruc :: !FrankCondon
               ,getCoeffCI :: ![[[Double]]]
                 } deriving Show
                                 
data CommonBlock = CB  {
                       getStates  :: Int,
                       getDt      :: Double,
                       getSubStep :: Int,
                       getETot    :: Double,
                       getDeco    :: Double,
                       getVp      :: [Double],
                       getVpp     :: [Double],
                       getVppp    :: [Double],
                       getp       :: [[Double]],
                       getpp      :: [[Double]],
                       getppp     :: [[Double]],
                       getRandom  :: [Double]
                       } deriving Show                 

-- | Data containing the state of the Thermostat
data Thermo = Thermo {
                     thQ1  :: !Double
                    ,thQ2  :: !Double
                    ,thVx1 :: !Double
                    ,thVx2 :: !Double
                 } deriving Show
                 
defaultMol :: Molecule 
defaultMol = Molecule zero zero zero zero [] [] [] (Left S0) (FC VU.empty V.empty) []
  where zero = R.fromListUnboxed  (Z:.(1::Int)) [0]
        
        
defaultCB :: CommonBlock 
defaultCB = CB 0 0 0 0 0 [] [] [] [[]] [[]] [[]] []
        
lookupLabel :: [(Label,a)] -> Label -> a
lookupLabel mapa l = 
   case DL.lookup l mapa of
        Nothing -> error $"The Label" DL.++ l DL.++ "was not found in the input file"
        Just xs -> xs        