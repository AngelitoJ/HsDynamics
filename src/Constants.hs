
module Constants where

-- =================> Constant <=================
ehKcal = 627.509         -- from Hartrees to Kcal/mol
kb = 3.1668152037e-6     -- Boltzmann constant in au hartree/Kelvin
amu = 1822.88851543      -- relation N/e
h = recip (2*pi)         -- planck's constant
a0 = 0.529177249         -- from bohr to angstrom
au_time = 2.41888432e-2  -- from femtoseconds to au
auN = 8.2387225e-8       -- atomic unit of force  a.u.
ev2kcal = 23.061         -- ElectronVolt to kcal/mol

mapCharge2Label :: [(Int,String)]
mapCharge2Label = [(1,"H"),(2,"He"),(3,"Li"),(4,"Be"),(5,"B"),
              (6,"C"),(7,"N"),(8,"O"),(9,"F"),(10,"Ne"),
              (11,"Na"),(12,"Mg")]

charge2Label :: Int -> String
charge2Label i = case lookup i mapCharge2Label of
                      Just l  -> l
                      Nothing -> error "Sorry, but I don't know what atoms are you talking about"