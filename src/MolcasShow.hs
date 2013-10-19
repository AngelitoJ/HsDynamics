-- ------------------------------------------------------------------------------------
--  Molcas 2010 (7.6) show instances for test porpouses
--  Authors: Angel Alvarez, Felipe Zapata
--
--  
-- ------------------------------------------------------------------------------------

module MolcasShow where

import qualified Data.List.Split as DLS
import Molcas

instance Show MolBlock where
    show (TextData s)             = "Text    : " ++ s
    show (LicenseData s)          = "License : " ++ s
    show (LogoData s)             = "Logo    : " ++ s
    show (CopyRightData s)        = "Copy Right : " ++ s
    show (InputData s)            = "Input   : " ++ s
    show (ProjectData p f s x)    = "Project info:\n\tName: " ++ p ++ "\n\tSubmitted: " ++ f ++ "\n\tScratch: " ++ s ++ "\n\tOuputs: " ++ x
    show (ModuleAuto n s modules) = "Module : " ++ n ++ " " ++ s ++ "\n" ++ concat (map show modules)

instance Show ModuleData where
    show (ModuleData n s items)   = "Module : " ++ n ++ " " ++ s ++ "\n" ++ concat (map show items)
    
instance Show SectionData where
    show (WhateverData s)                   = "Module Data: " ++ (show s) ++ "\n"
    show (RasSCFData s)                     = "RasSCF Data: " ++ (show s) ++ "\n"
    show (RasSCFCI r e ci)                  = "RasSCF CI-Coefficients for root: " ++ (show r) ++ " Energy: " ++ (show e) ++ " info:" ++ (show ci) ++ "\n"
    show (RasSCFRE roots)                   = "RasSCF Final Energies:" ++ (show roots)++"\n"
    show (AlaskaData s)                     = "Alaska Data: " ++ (show s) ++ "\n"
    show (AlaskaMolecularGradients r grads) = "Alaska Molecular Gradients [" ++ r ++ "]:" ++ (showMolGrads grads) ++ "\n"


-- Make a Alaska Molecular Gradients listing suitable for humans...
showMolGrads :: [(String,Char,Double)] -> String
showMolGrads list = concat $ map showGrad $ chunkedlist
    where
        chunkedlist = DLS.chunksOf 3 list
        showGrad [(a1,x,vx),(a2,y,vy),(a3,z,vz)] = "\n\tAtom " ++ a1 
            ++ "  " ++ [x] ++ " = " ++ (show vx)
            ++ ", " ++ [y] ++ " = " ++ (show vy)
            ++ ", " ++ [z] ++ " = " ++ (show vz)
