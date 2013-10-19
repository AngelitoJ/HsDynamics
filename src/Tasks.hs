-- Task management for parsing jobs
-- @2013 Angel Alvarez

-- {-# LANGUAGE FlexibleInstances #-}

module Tasks where

import Data.Function
import Data.List

-- Processor Affinity type allows to map a list of task to certain processor cores, Affinities entities only 
-- order and compare on the processor fields so "Ord" and "Eq" instances need to to be redefined to allow this 
-- behavour
--
-- Using GroupBy (compare `on` fst) $ SortBy (compare `on` fst) can achieve similar results without the need 
-- to use newtypes besides the lack of elegance.

-- A similiar mapping can be obtained mapping tasks modulo capabilities


newtype Affinity a b = Affinity (a, b) 
--   deriving (Show)

instance (Show a, Show b) => Show (Affinity a b) where
    show (Affinity (p, t)) = "Affinity " ++ show (p, t)
    
instance (Eq a, Num a, Eq b) => Eq (Affinity a b) where
    (==) (Affinity (p0, _)) (Affinity (p1, _)) = p0 == p1          -- We only check for processor element

instance (Ord a, Num a, Ord b) => Ord (Affinity a b) where
    compare (Affinity (p0,_)) (Affinity (p1,_)) = compare p0 p1    -- We only check for processor element

-- Public Constructor to avoid users messing up with internal details 
toAff :: (a, b) -> (Affinity a b)
toAff = Affinity

-- Public Deconstructor
fromAff :: (Affinity a b) -> b
fromAff (Affinity (a, b)) = b


-- Map files to cores, and make sublists so every core gets a file list to compute
-- mapfilestocores :: Int  -> [ String ]-> [[ String ]]
-- mapfilestocores cores files = map (\x -> map fromAff x) $ group $ sort $ map toAff $ zip (cycle corelist ) files
--     where
--         corelist = [0..cores-1]

-- Map files to cores, and make sublists so every core gets a file list to compute
-- f2c :: Int  -> [ String ]-> IO ( [[ String ]] )
-- f2c cores files = do
--     map (\x -> map fromAff x) $ group $ sort $ map toAff $ zip (cycle corelist ) files
--     where
--         corelist = [0..cores-1]



-- Test functions for new type and basic instances

-- Full Polymorphic test function showing the overloaded instances at full pace
-- test :: (Enum a, Eq a, Num a) => [[Affinity a [Char]]]
-- test = group $ map Affinity $ zip (cycle [0..3]) ["a","b","c","d","f"]

-- Usual groupping using just overloading when needed
-- test2 :: [[Affinity Int [Char]]]
-- test2 = group $ map Affinity $ sort $ zip (cycle [0..3]) ["a","b","c","d","f"]

-- Semi complete groupping test using constructors and destructor to change common behaviors
-- test3 :: [[ String ]]
-- test3 = map (\x -> map fromAff x) $ group $ map toAff $ sort $ zip (cycle [0..3]) ["a","b","c","d","f"]

-- IO task designed to run a long lasting script, we dont care about return codes or such

-- processOne :: IO ()
-- processOne = do
--     putStrLn "Starting busy wait script\n"
--     threadDelay 2000000
--     _ <- system "./busy.sh"                     -- 
--     putStrLn "Finished busy wait script\n"
--     return ()

-- IO task carring out main parsec processing
-- processTwo :: Options -> [FilePath] -> IO ()
-- processTwo o f = do
--     putStrLn "Starting parsing task\n"
--     mapM_ (processFile o) f

--  A simple IO processs manager, manages Asyncs to acomplish several tasks...
-- processZero ::  Options -> [FilePath] -> IO ()
-- processZero o f = do
--     putStrLn "Starting background tasks\n"
--     concurrently (processOne) (processTwo o f)
--     putStrLn "Finished al tasks\n"
    
