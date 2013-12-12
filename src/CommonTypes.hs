{-# Language DeriveFunctor,FlexibleInstances , TemplateHaskell #-}


module CommonTypes where

import Control.Lens
import Data.Array.Repa as R hiding ((++))
import Data.ByteString.Char8 (ByteString(..))
import Data.Complex
import Data.Functor.Identity
import Data.List as DL
import qualified Data.List.Split as DLS
import qualified Data.Map as M
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector as V
import Text.Parsec
import Text.Parsec.ByteString
import Text.Printf


-- =============================================================================
data Singlet = S0 | S1 | S2 | S3 | S4 deriving (Show,Enum,Eq,Read)

data Triplet = T1 | T2 | T3 | T4 | T5 deriving (Show,Enum,Eq,Read)

data TheoryLevel = CASSCF (Int,Int) Int String  | HF | Unspecified   -- (number Electrons, Orbital) relaxRoot further commands

instance Show TheoryLevel where
    show (CASSCF (electrons,orbitals) rlxroot s) = "CASSCF" ++ "(" ++ show electrons ++ "," ++ show orbitals ++ ",NRoot=" ++ show rlxroot ++ s ++ ")" 
    show HF                                      = "HF"
    show Unspecified                             = "Unspecified"
 
data MolcasInput a =  Command  a |  Gateway a | Seward a | RasSCF Int a a | MCLR a | Alaska a deriving (Functor)

instance Show (MolcasInput String) where 
  show (Command x)    = ">> "        ++ x ++ "\n"
  show (Gateway x)    = "&Gateway\n" ++ x
  show (Seward  x)    = "&Seward\n"  ++ x
  show (RasSCF n x y) = "&Rasscf\n"  ++ x ++ "rlxroot=" ++ show n ++ y
  show (MCLR    x)    = "&Mclr\n"   ++ x
  show (Alaska  x)    = "&Alaska\n"  ++ x
 
data Job = Gaussian (TheoryLevel,Basis) | Interpolation | Molcas [MolcasInput String] | MolcasTinker [(Label,Int)] Command 
          | Palmeiro Connections [FilePath]|Quadratic | HaskellAbInitio deriving Show 

-- Internal Coordinates types 
data  InternalCoord = Bond !Int !Int | Angle !Int !Int !Int | Dihedral !Int !Int ! Int !Int deriving Show                          
type Connections = V.Vector InternalCoord
type Force     = Array U DIM1 Double 
type Grad      = Array U DIM1 Double
type Hess      = Array U DIM2 Double
type Internals = VU.Vector Double  

type Anchor      = [Int]
type Args        = String
type Basis       = String
type Command     = String
type Coordinates = Array U DIM1 Double
type DT          = Double -- step of integration
type Energy      = Double
type Label       = String
type Lines       = String
type Masses      = Array U DIM1 Double
type MatrixCmplx = [Complex Double]
type Name        = String
type Parameters  = String
type Project     = String
type ReactionCoordinate = M.Map Int (Int,Internals)
type Roots       = Int    -- Number of roots of the system
type Temperature = Double
type TotalEnergy = Double
type Time        = Double
type VelocityXYZ = [Double]
type XYZ         = [Double]

data InitialDynamics = InitialDynamics {
                                        _getInitialState :: Singlet
                                       ,_getTime         :: Double
                                       ,_getdt           :: Double
                                       ,_getTheory       :: TheoryLevel
                                       ,_getBasis        :: Basis
                                       ,_getExtForceMod  :: Double 
                                       ,_getForceAnchor  :: [Int]
                                       ,_getProject      :: String
                                       }    deriving Show

makeLenses ''InitialDynamics                                        
                                       
-- ==================> Internal Coordinates data types <============
data EnergyDerivatives = EnergyDerivatives !Energy !Grad !Hess deriving Show

data FrankCondon = FC !Internals !Connections deriving Show



-- ========================> Types for Parsing <======================================
data BlockType = GauBlock | MolBlock deriving Show     

data EigenBLock = EigenBLock { getRootEnergy :: !Double 
                              ,getRootCoeff ::  ![Double] }  deriving Show

data GauBlock = 
      IGauBlock Label Int [Integer]     -- A GauBlock containing one or more integer values
    | RGauBlock Label Int [Double]      -- A GauBlock containing one or more real values
    | TGauBlock String                  -- A GauBlock containing just unformatted text

    
data GaussLog = GaussLog [EigenBLock] [[Double]] deriving Show

                                       
-- ==========> Molcas and Gaussian parser Types <=========                                      
type BlockParser   = MyParser MolState MolBlock
type SectionParser = MyParser MolState SectionData 

-- A parsecT based parser carring user supplied state of type "u" and returning data of type "a"
type MyParser u a = ParsecT ByteString u Identity a


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
    | AlaskaMolecularGradients String [(String,Char,Double)] [[Double]] -- Alaska Molecular Gradients and optional Gradient after ESPF

-- Common state record to all Molcas parsers
data MolState = MolState 
    { 
          blockLevel   :: [BlockParser]                       -- Block level Molcas blocks
        , sectionLevel :: [SectionParser]                     -- Section Level Molcas Blocks
    }

         
-- =================> Quantum Mechanics and Molecular Mechanics Atoms <=====
data AtomMM = AtomMM {_numberMM :: Int, _labelMM :: Label , _xyzMM :: XYZ, _parametersMM :: Parameters} deriving Show

makeLenses ''AtomMM 

data AtomQM = AtomQM Label XYZ deriving Show

data AtomXYZ = Atom Label XYZ VelocityXYZ deriving Show


-- =================> Types for Dynamics <===============

data Molecule = Molecule {
                _getCoord   :: !(R.Array R.U R.DIM1 Double)
               ,_getVel     :: !(R.Array R.U R.DIM1 Double)
               ,_getForce   :: !(R.Array R.U R.DIM1 Double)
               ,_getMass    :: !(R.Array R.U R.DIM1 Double)
               ,_getAtoms   :: ![String]
               ,_getEnergy  :: ![[Double]]
               ,_getDervEn  :: ![EnergyDerivatives]
               ,_getElecSt  ::  Either Singlet Triplet
               ,_getFCStruc :: !FrankCondon
               ,_getCoeffCI :: ![[[Double]]]
                 } deriving Show
                 
makeLenses ''Molecule                 
                 
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
                 

-- ================> Instances Molcas <=================
instance Show MolBlock where
    show (TextData s)             = "Text    : " ++ s
    show (LicenseData s)          = "License : " ++ s
    show (LogoData s)             = "Logo    : " ++ s
    show (CopyRightData s)        = "Copy Right : " ++ s
    show (InputData s)            = "Input   : " ++ s
    show (ProjectData p f s x)    = "Project info:\n\tName: " ++ p ++ "\n\tSubmitted: " ++ f ++ "\n\tScratch: " ++ s ++ "\n\tOuputs: " ++ x
    show (ModuleAuto n s modules) = "Module : " ++ n ++ " " ++ s ++ "\n" ++ concatMap show modules

instance Show ModuleData where
    show (ModuleData n s items)   = "Module : " ++ n ++ " " ++ s ++ "\n" ++ concatMap show items

instance Show SectionData where
    show (WhateverData s)                   = "Module Data: " ++ (show s) ++ "\n"
    show (RasSCFData s)                     = "RasSCF Data: " ++ (show s) ++ "\n"
    show (RasSCFCI r e ci)                  = "RasSCF CI-Coefficients for root: " ++ (show r) ++ " Energy: " ++ (show e) ++ " info:" ++ (show ci) ++ "\n"
    show (RasSCFRE roots)                   = "RasSCF Final Energies:" ++ (show roots)++"\n"
    show (AlaskaMolecularGradients r grads gradESPF)
                                            = "Alaska Molecular Gradients [" ++ r ++ "]:" ++ (showMolGrads grads) ++ "\n" ++ 
                                              (if null gradESPF then []  else "Gradient after ESPF\n" ++ show gradESPF)
    show (AlaskaData s)                     = "Alaska Data: " ++ (show s) ++ "\n"
                                                                
                                                         

    
-- Make a Alaska Molecular Gradients listing suitable for humans...
showMolGrads :: [(String,Char,Double)] -> String
showMolGrads list = concatMap showGrad $ chunkedlist
    where
        chunkedlist = DLS.chunksOf 3 list
        showGrad [(a1,x,vx),(a2,y,vy),(a3,z,vz)] = "\n\tAtom " ++ a1
            ++ "  " ++ [x] ++ " = " ++ (show vx)
            ++ ", " ++ [y] ++ " = " ++ (show vy)
            ++ ", " ++ [z] ++ " = " ++ (show vz)


defaultMol :: Molecule 
defaultMol = Molecule zero zero zero zero [] [] [] (Left S0) (FC VU.empty V.empty) []
  where zero = R.fromListUnboxed  (Z:.(1::Int)) [0]
        
        
defaultCB :: CommonBlock 
defaultCB = CB 0 0 0 0 0 [] [] [] [[]] [[]] [[]] []
        
lookupLabel :: [(Label,a)] -> Label -> a
lookupLabel mapa l = 
   case DL.lookup l mapa of
        Nothing -> error $"The Label  " ++ l ++ "  was not found in the input file"
        Just xs -> xs        