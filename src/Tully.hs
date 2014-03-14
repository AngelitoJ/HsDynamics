

module Tully (
              initialAMTX
             ,tullyHS
             ) where

import Control.Applicative             
import Control.Arrow ((&&&))   
import Control.Concurrent (forkIO)
import Control.Exception (assert)
import Control.Lens  hiding (Index)
import Data.List as L
import Data.Array.Repa as R hiding ((++))
import Data.Complex
import Data.List.Split
import Data.Maybe
import qualified System.Random as SR
import Data.Function (on)
import Text.Printf


-- ===> Internal Modules <==
import CommonTypes
import Dynamics

-- ===============> Data Types <=================
type Coefficients = [[Double]]
type SubStepIndex = Int
type MatrixReal = [Double]
type Index = Int
type RealValue = Double
type Vvector = [Double]
type RlxRoot = Int


-- =============> Testing Values <================
{-
a=[-0.802440331997109    , -0.575868486355278      , -4.727377639449486E-003 , -1.061564961986605E-002 ,  6.174697331292390E-002 ,  5.242139273240287E-002 ,  3.057781545313075E-003 ,  1.554930566599939E-002 , -3.062390943705189E-003 , -7.960272362267257E-002 , -1.845966029448955E-003 ,  1.334119944488831E-002 ,  2.812495057039225E-004 ,  5.881227869310120E-002 , -7.477199881204565E-002 , -1.394767489469036E-002 , -8.668648147390514E-004 ,  3.248107205578303E-004 ,  4.027173291624844E-002 , -1.112423889766913E-002 ,  0.567648542467870      , -0.792058640754844      , -2.024363905458087E-003 , -0.105034595125135      , -0.167507441152981      ,  2.677403176120711E-002 ,  1.081238429166823E-003 , -1.391804620903221E-002 , -2.326891137513524E-003 , -5.933558139209187E-002 ,  3.186083828273563E-004 , -1.571171294301005E-002 , -8.436117974534084E-003 , -3.974417352604644E-002 , -3.459530783059087E-002 , -3.149404342856874E-002 , -2.615191944975547E-003 , -1.880121868571484E-003 ,  5.163218519090958E-002 , -1.
107184220949405E-002 ,  5.080787484372048E-002 , -8.955367014682171E-002 ,  5.404836165526769E-003 ,  0.990489102592445      , -2.247416535264812E-002 ,  6.700682151808082E-003 ,  3.101108582367395E-002 , -1.988345486615977E-003 ,  5.169656199141524E-002 , -5.683226900788068E-003 , -2.081939176788467E-003 , -2.337789415100215E-004 ,  4.763876870337029E-002 , -2.745413330988324E-003 , -4.693871642877469E-003 , -1.636533271701089E-003 ,  1.947063841969837E-002 ,  3.705513279264840E-002 ,  4.569423971948600E-003 , -1.098167594643677E-003]-- {{{
b=[-0.825300493917011      ,-0.543879261908658      ,-7.482993749622175E-003 ,-7.098395645823079E-003 , 5.673057014923283E-002 , 5.716416938811585E-002 , 5.000828353357999E-003 , 1.520066688885394E-002 ,-3.733489320256863E-003 ,-7.581468064748702E-002 ,-3.739033221426126E-003 , 1.317755764159942E-002 , 2.088795550784969E-003 , 5.696902783331610E-002 ,-7.340815433591982E-002 ,-1.383491534062133E-002 ,-1.007898827172710E-003 , 2.602761821374385E-004 , 3.636393324564140E-002 ,-1.109255570174482E-002 , 0.532732309241082      ,-0.808554269231076      ,-1.469450139485942E-003 ,-0.155465136807209      ,-0.166743894822927      , 2.598442318235850E-002 , 1.646122755480915E-003 ,-1.235477917135902E-002 ,-4.553987528458254E-003 ,-5.574117970150035E-002 , 1.638756729309313E-003 ,-1.446447128261960E-002 ,-9.101007158706733E-003 ,-3.428061783231465E-002 ,-3.322375490755193E-002 ,-3.614026003392144E-002 ,-4.085593375395304E-003 ,-3.662134222010870E-003 , 4.963182293965447E-002 ,-1.138307783318109E-002 , 7.774743138333563E-
002 ,-0.130043642271938      , 5.712036173459300E-003 , 0.983846968246718      ,-3.247021357500395E-002 , 1.555290007983219E-002 , 2.787996534720727E-002 ,-2.807168513952760E-003 , 5.371075629525862E-002 ,-1.192539985159257E-002 ,-2.036664731115888E-003 ,-8.543778752505781E-004 , 4.519432828944159E-002 ,-3.945950140286140E-003 ,-7.256587540268542E-003 ,-9.335332134310540E-004 , 1.803971817221128E-002 , 3.840340984937011E-002 , 5.327039044334785E-003 ,-1.710873196627125E-003]
c=[-0.843449984658543     ,-0.516652157938039     ,-1.089492957822513E-002,-8.696132428339920E-003, 4.981652698949847E-002, 6.205213163354716E-002, 8.650926017709035E-003, 1.489567129117041E-002,-3.525182903182093E-003,-7.215693969432212E-002,-6.362979566917558E-003, 1.347139206612942E-002, 2.918353457808912E-003, 5.329183501896391E-002,-7.208745155403193E-002,-1.316296909176547E-002,-1.450640002840407E-003,-6.166156924805577E-005, 3.261919148912029E-002,-1.067419518537761E-002, 0.499045904031500     ,-0.811046081695021     , 9.737858684978825E-004,-0.238938236238964     ,-0.163025311798853     , 2.473867207459687E-002, 2.704570897295054E-003,-1.102212004963909E-002,-8.813570008882898E-003,-5.136680237439385E-002, 3.640717893187751E-003,-1.338487169512298E-002,-9.112285129240262E-003,-2.782634977216089E-002,-3.292128054070918E-002,-3.924474144136334E-002,-6.545548471047820E-003,-7.135067586666743E-003, 4.629130583965130E-002,-1.097172220212339E-002, 0.116027790150284     ,-0.202351040639085     , 7.
549231016298813E-003, 0.966803487793439     ,-4.561323412114333E-002, 2.845123472852012E-002, 2.566357151982683E-002,-3.828148174510277E-003, 5.833014511393888E-002,-2.055843606840661E-002,-1.560360366319303E-003,-1.238210222827141E-003, 3.785537348990355E-002,-5.572738008838105E-003,-1.151329664538233E-002,-5.375364080672151E-004, 1.658199360225748E-002, 3.880817833476892E-002, 6.936396484634204E-003,-2.686331913097163E-003]-- }}}
states=3 :: Int
dt   = 20.65
nSubStep = 200 :: Int
eTot = -93.5305237411652
deco = 0.1
vp=[-0.938383638E+02,-0.936493854E+02,-0.934944314E+02]
vpp=[-0.938444231E+02,-0.936306617E+02,-0.934885032E+02]
vppp=[-0.938443372E+02,-0.936089491E+02,-0.934829116E+02]
p=chunksOf  (floor $ fromIntegral (length a) / fromIntegral states) a
pp=chunksOf  (floor $ fromIntegral (length a) / fromIntegral states) b
ppp=chunksOf  (floor $ fromIntegral (length a) / fromIntegral states) c
aMatrix=[3.631964048315513E-005 :+ 1.091656294664714E-007 , 1.690700266553072E-003 :+ 6.190927351558496E-003 , 1.752692661446557E-004 :+ (-2.280270230623447E-004) , 1.745534742823029E-003 :+ (-6.326677447060865E-003) , 0.998803839851981 :+ 0.000000000000000E+000 , (-1.947873148665928E-002) :+ (-3.092162117050096E-002) , 1.790857066777628E-004 :+ 2.339424354093770E-004 , (-1.951268569833850E-002) :+ 3.098670398771564E-002 , 1.159840507535778E-003 :+ (-2.021760112245944E-007)]
rlxRoot=2 :: Int-}


-- ======> <===========
tullyHS :: DT -> MatrixCmplx -> Int -> Molecule -> IO (Molecule, MatrixCmplx)
tullyHS  dt aMatrix step mol = do
    let energies      = mol ^. getEnergy
    if length energies <3 then return (mol,aMatrix)
                          else tullyAction energies
                          
  where tullyAction energies = do                          
           lo <- SR.randomIO
           let eTot          = calcTotalEnergy mol
               seed          = SR.mkStdGen lo
               random        = take nSubStep $ filter (> 10**(-4)) $ aleGenRan seed                
               coeffs        = mol ^. getCoeffCI 
               [p,pp,ppp]    = coeffs
               [vp,vpp,vppp] = energies
               states        = 2
               nSubStep      = 200
               deco          = 0.1
               correctedp    = correctSignsByLuisma p pp vp rlxRoot
               rlxRoot       = succ $ calcElectSt mol 
               cb            = CB states dt nSubStep eTot deco vp vpp vppp correctedp pp ppp random
               f1 = "Tully Haskell Routine for Dynamic step:" ++ (show $ succ step) ++ "\n"
               f2 = "Substeps: " ++ (show nSubStep) ++ "\n"
               f3 = "Decoherence Factor: " ++ (show deco) ++ "\n"
               f4 = "Step delta: "  ++ (show dt) ++ "\n"
           appendFile "TullyOutput" $ L.foldl1' (++) [f1,f2,f3,f4]
           (newRelax,newAmatrix) <- driver aMatrix rlxRoot 1 cb mol
           if newRelax == rlxRoot then return (mol,newAmatrix) 
                                  else  do
                                        let intSt = pred newRelax
                                            newSt    =  toEnum  $ intSt -- Tully Modules number S0 with 1 Dynamics with 0 
                                            newCoeff = correctedp : (tail coeffs)
                                            currentEnergy = let xs = mol^.getEnergy . to head   in xs !! intSt
                                            newVel = scaleVelocity eTot currentEnergy mol
                                            newMol = mol & getCoeffCI .~ newCoeff 
                                                         & getElecSt  .~ (Left newSt) 
                                                         & getVel     .~ newVel
                                        appendFile "result.out" $ "Hop to Root: " ++ (show newSt) ++ "\n" 
                                        return (newMol,newAmatrix)
                                     
initialAMTX :: Molecule -> MatrixCmplx
initialAMTX mol = [fun i j | i<- [0..st], j <- [0..st] ]
  where st = calcElectSt mol
        fun x y  = if all (==st) [x,y] then 1.0 :+ 0.0
                                       else 0.0 :+ 0.0
                                         
                                     
scaleVelocity :: Double -> Double -> Molecule -> Array U DIM1 Double
scaleVelocity etot ep mol = let [vs,ms] = fmap (mol^.) [getVel,getMass]
                                ek = calcEk vs ms
                                sub = abs $ etot - ep   
                                s = sqrt (sub/ek)
                            in computeUnboxedS $ R.map (*s) $  vs                               
 
driver :: MatrixCmplx -> Int -> Int -> CommonBlock  -> Molecule -> IO (Int,MatrixCmplx)
driver aMat rlxRt step cb mol | step <= getSubStep cb = do
        let
           vMatS   = matrixVSub cb step
           dMatS   = matrixDsub cb step
           aDT     = calculateAdt cb vMatS dMatS aMat
           aMat'   = integrateA cb aMat aDT
--           probab  = probability cb aMat' dMatS rlxRt 
        probab  <- probabilityIO cb aMat' dMatS vMatS rlxRt
        let
           newRoot = checkHop cb probab rlxRt step mol        
           states  = getStates cb
           aMat''  = perGranCorr cb aMat' vMatS newRoot   
           f1      = "\nStep: " ++ (show step) ++ "\n"
           f2      = "Matrix D:\n"
           f3      = printWell dMatS states
           f4      = "Matrix V:\n"
           f5      = printWell vMatS states
           f6      = "Matrix Adt:\n"
           f7      = printWellC aDT states
           f8      = "Random: " ++ (show $ (getRandom cb) !! (pred step)) ++ "\n"
           f9      = "Probabilities: "
           f10     = printWell probab states
           f96     = "Hop after checking Root: " ++ show (newRoot) ++ "\n"
           f98     = "Matrix A after Integration:\n"
           f99     = printWellC aMat' states
           f11     = "Matrix A after Persico/Granucci correction:\n"
           f12     = printWellC aMat'' states              
           string  =  L.foldl1' (++) [f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f96,f98,f99,f11,f12]       
        appendFile "TullyOutput" string       
        driver aMat'' newRoot (succ step) cb mol
      | otherwise = return (rlxRt, aMat)

correctSignsByLuisma :: Coefficients -> Coefficients -> Vvector -> RlxRoot -> Coefficients
correctSignsByLuisma p pp vp rlxrt = fmap (correctOneByLuisma p pp) $ goodCouplesByLuisma p pp $ rightOrderByLuisma vp rlxrt

correctOneByLuisma :: Coefficients -> Coefficients -> (Int, Int) -> [Double]
correctOneByLuisma p pp (a,b) = let 
                sign = signum $ scalarProd (p!!a) (pp!!b)
                in fmap (fmap (*sign)) p!!a

rightOrderByLuisma :: Vvector -> RlxRoot -> [Int]
rightOrderByLuisma vp rlxrt = let
    tupla            = zip3 [0..] vp (fmap abs (fmap (+(-(vp!!(rlxrt-1)))) vp))
    orderedTupla     = sortBy (compare `on` trd') tupla
    in fmap fst' orderedTupla

goodCouplesByLuisma :: Coefficients -> Coefficients -> [Int] -> [(Int, Int)]
goodCouplesByLuisma p pp rightOrder = sortBy (compare `on` fst) $ snd $ foldl go acco rightOrder
     where tupl1 = zip p [0..]
           tupl2 = zip pp [0..]
           acco  = (tupl2,[])
           go (pps,resTupl) ip = let
               psVal  = fmap fst tupl1 
               vi     = psVal!!ip
               prosca = fmap (scalarProdT vi) pps
               ipp    = snd $ maximum prosca 
               newxs  = filter (\(x,ip) -> ip /= ipp)  pps
               in (newxs,(ip,ipp):resTupl)

scalarProdT :: [Double] -> ([Double],Int) -> (Double,Int)
scalarProdT a b = (abs $ sum $ L.zipWith (*) a (fst b), (snd b))

scalarProd :: [Double] -> [Double] -> Double
scalarProd a b = sum $ L.zipWith (*) a b

fst' (a,_,_) = a
snd' (_,b,_) = b
trd' (_,_,c) = c

createA :: CommonBlock -> RlxRoot -> MatrixCmplx
createA cb rlxRoot = let
      states     = getStates cb
      zeroesA    = replicate states $ replicate states (0 :+ 0)
      rightIndex = states*(rlxRoot-1) + rlxRoot
      matrixA    = (take (rightIndex-1) $ concat zeroesA) ++ [1 :+ 0] ++ (drop rightIndex $ concat zeroesA) 
      in matrixA

subAtIndex :: [Index] -> RealValue -> MatrixReal -> MatrixReal
subAtIndex is n xs = fmap (\(i,x)-> if elem i is then n else x) $ zip [0..] xs

subAtIndex' :: [Index] -> [RealValue] -> MatrixReal -> MatrixReal
subAtIndex' is doubles xs = let
             rightIndex x = maybe 0 id $ elemIndex x is 
             in fmap (\(i,x)-> if elem i is then doubles!!(rightIndex i) else x) $ zip [0..] xs

createD :: CommonBlock -> MatrixReal
createD cb = let
        p        = getp cb
        pp       = getpp cb
        states   = getStates cb
        dt       = getDt cb
        elementD pE ppE = -(sum $ L.zipWith (*) pE ppE )/dt
        currentI = concatMap (replicate states) p
        previouI = concat $ replicate states pp
        matrixD' = L.zipWith elementD currentI previouI
        diagonalI= L.zipWith (+) [0..states] (take states (iterate (+states) 0))
        matrixD  = subAtIndex (diagonalI) 0.0  matrixD'
        in matrixD

creatDhalf :: CommonBlock -> Coefficients -> Coefficients -> MatrixReal
creatDhalf cb p pp = let
        dt       = getDt cb
        states   = getStates cb
        scalPr x y = sum $ L.zipWith (*) x y
        element32 a b c d = ((scalPr a b)-(scalPr c d))/(2 * dt)
        aaabbbccc = concatMap (replicate states) 
        abcabcabc = concat . replicate states
        currentI  = aaabbbccc pp -- dddeeefff
        previouI  = abcabcabc p  -- abcabcabc
        currentII = aaabbbccc p  -- aaabbbccc
        previouII = abcabcabc pp -- defdefdef
        matrixD32 = L.zipWith4 element32 currentI previouI currentII previouII
        in matrixD32

creatDhalf' :: CommonBlock -> Coefficients -> Coefficients -> MatrixReal
creatDhalf' cb p pp = let
        dt       = getDt cb
        scalPr w z = sum $ L.zipWith (*) w z
        xs =  L.zipWith (++) [[x,y] | x<-pp, y<-p] [[x,y] | x<-p, y<-pp]
        in fmap (\[a,b,c,d] -> ((scalPr a b)-(scalPr c d))/(2 * dt)) xs

matrixDExtrSlope :: CommonBlock -> MatrixReal
matrixDExtrSlope cb = let
         p       = getp cb
         pp      = getpp cb
         ppp     = getppp cb
         dt      = getDt cb 
         d32     = creatDhalf cb pp ppp 
         d12     = creatDhalf cb p pp 
         in L.zipWith (\x y -> (x-y)/dt) d12 d32

matrixDExtrInter :: CommonBlock -> MatrixReal
matrixDExtrInter cb = let
         p       = getp cb
         pp      = getpp cb
         in creatDhalf  cb p pp

matrixDsub :: CommonBlock -> SubStepIndex -> MatrixReal
matrixDsub cb n = let
        p       = getp cb
        pp      = getpp cb
        ppp     = getppp cb
        dt      = getDt cb
        nSubStep = getSubStep cb 
        slope   = matrixDExtrSlope cb
        inter   = matrixDExtrInter cb
        counter = ((fromIntegral n-1) * dt)/(fromIntegral nSubStep) - dt/2
        in L.zipWith (+) inter (fmap (* counter) slope)

matrixVSub :: CommonBlock -> SubStepIndex -> MatrixReal
matrixVSub cb n = let
        vp       = getVp cb
        vpp      = getVpp cb
        dt       = getDt cb
        states   = getStates cb
        nSubStep = getSubStep cb
        slope    = L.zipWith (\x y -> (x-y)/dt) vp vpp
        inter    = vpp
        counter  = ((fromIntegral n-1) * dt)/(fromIntegral nSubStep)
        diagonal = L.zipWith (+) inter (fmap (* counter) slope)
        zeroeV   = replicate (states*states) (0.0)
        rightInd = L.zipWith (+) [0..states] (take states (iterate (+states) 0))
        in subAtIndex' rightInd diagonal zeroeV

calculateAdt :: CommonBlock -> MatrixReal -> MatrixReal -> MatrixCmplx -> MatrixCmplx
calculateAdt cb vMat dMat aMat = let
        states   = getStates cb
        columnI n = fmap (+n) $ take states $ iterate (+states) 0
        rowI n = fmap (+(n*states)) [0..states-1]
        rightTupla x = (div x states, mod x states)
        oneILtwoLJ one two (i,j) = sum $ L.zipWith (*) (fmap (one !!) $ rowI i) (fmap (two !!) $ columnI j)
        oneLJrwoIL one two (i,j) = sum $ L.zipWith (*) (fmap (one !!) $ columnI j) (fmap (two !!) $ rowI i) 
        dMatC = fmap (\x -> x :+ 0.0) dMat
        vMatC = fmap (\x -> x :+ 0.0) vMat
        first = fmap (\x -> oneILtwoLJ aMat vMatC (rightTupla x)) [0..(states*states-1)]
        second= fmap (\x -> oneLJrwoIL aMat vMatC (rightTupla x)) [0..(states*states-1)]
        third = fmap (\x -> oneILtwoLJ aMat dMatC (rightTupla x)) [0..(states*states-1)]
        fourth= fmap (\x -> oneLJrwoIL aMat dMatC (rightTupla x)) [0..(states*states-1)]
        ii    = 0.0 :+ 1.0
        in L.zipWith4 (\a b c d -> (ii * a) - (ii * b) + c - d) first second third fourth

integrateA :: CommonBlock -> MatrixCmplx -> MatrixCmplx -> MatrixCmplx
integrateA cb aMat aMatDT = let
              dt       = getDt cb
              nSubStep = getSubStep cb
              rightDT  = dt / (fromIntegral nSubStep) 
              in L.zipWith (\x y -> (x) + ((y)*(rightDT :+ 0.0))) aMat aMatDT

probability :: CommonBlock -> MatrixCmplx -> MatrixReal -> MatrixReal -> RlxRoot -> MatrixReal
probability cb aMat dMat vMat rlxRt = let 
    states           = getStates cb
    dt               = getDt cb
    nSubStep         = getSubStep cb
    fromArrtoMat v   = [if i==j then (v!!i :+ 0.0) else (0.0 :+ 0.0) | i <- [0.. pred (length v)], j <- [0..pred (length v)]] 
    dMatC            = fmap (\x -> x :+ 0.0) dMat
    vMatC            = fmap (\x -> x :+ 0.0) vMat
    matrixB          = L.zipWith3 (\x y z -> 2 * (imagPart ((conjugate $ x) * z)) - (2 * (realPart $ (conjugate $ x) * y ))) aMat dMatC vMatC
    rightIndex (n,m) = n*states + m
    columnI n        = fmap (+n) $ take states $ iterate (+states) 0
    bColumn          = fmap (matrixB !!) $ columnI (rlxRt-1)
    aTemTemp         = realPart $ aMat !! (rightIndex (rlxRt-1,rlxRt-1))
    rightDT          = dt / (fromIntegral nSubStep)
    probab'          = fmap (\x -> (x*rightDT)/aTemTemp) bColumn
    isPositive  x    = if x<0 then 0.0 else x
    probab''         = fmap isPositive probab'
    in makeRlxRtZero probab'' rlxRt
    

makeRlxRtZero :: [Double] -> Int -> [Double]
makeRlxRtZero xs i = let tuplas         = zip xs [0..]
                         ind            = pred i
                         tuplaFiltered  = fmap (\x -> if snd x /= ind then x else (0.0,ind)) tuplas
                     in fmap fst tuplaFiltered


probabilityIO :: CommonBlock -> MatrixCmplx -> MatrixReal -> MatrixReal -> RlxRoot -> IO(MatrixReal)
probabilityIO cb aMat dMat vMat rlxRt = do
    let states           = getStates cb
        dt               = getDt cb
        nSubStep         = getSubStep cb
        fromArrtoMat v   = [if i==j then (v!!i :+ 0.0) else (0.0 :+ 0.0) | i <- [0.. pred (length v)], j <- [0..pred (length v)]] 
        dMatC            = fmap (\x -> x :+ 0.0) dMat
        vMatC            = fmap (\x -> x :+ 0.0) vMat
        matrixB          = L.zipWith3 (\x y z -> 2 * (imagPart ((conjugate $ x) * z)) - (2 * (realPart $ (conjugate $ x) * y ))) aMat dMatC vMatC
        rightIndex (n,m) = n*states + m 
        columnI n        = fmap (+n) $ take states $ iterate (+states) 0
        bColumn          = fmap (matrixB !!) $ columnI (rlxRt-1)
        aTemTemp         = realPart $ aMat !! (rightIndex (rlxRt-1,rlxRt-1))
        rightDT          = dt / (fromIntegral nSubStep)
        probab'          = fmap (\x -> (x*rightDT)/aTemTemp) bColumn
        isPositive x     = if x<0 then 0.0 else x
        probab''         = fmap isPositive probab'
        probab'''        = makeRlxRtZero probab'' rlxRt
        f1 =  "A:"
        f2 = show aMat
        f3 = "B:"
        f4 = show matrixB
        f5 = "Values V"
        f6 = show vMat
        f7 = "dmat complex"
        f8 = show dMatC
        string  =  L.foldl1' (++) [f1,f2,f3,f4,f5,f6,f7,f8]
    appendFile "TullyOutput" string
    return $ probab'''


checkHop :: CommonBlock -> MatrixReal -> SubStepIndex -> RlxRoot -> Molecule -> RlxRoot
checkHop cb probab rlxRt step mol = let
         statesI   = pred (getStates cb)
         random    = assert (pred step >= 0) $ (getRandom cb) !! (step -1)
         sumProb   = fmap (\x -> sum $ take (x+1) probab) $ [0..statesI]
         etot      = calcTotalEnergy mol
         jumpTo    = dropWhile ((<random).fst) $ zip sumProb [1..]
         newRoot   = snd . head $ jumpTo
         energies  = mol ^. getEnergy . to head 
         jumpEnergy = energies !! (pred newRoot)
         in if jumpTo == [] then rlxRt else if jumpEnergy > etot then rlxRt else newRoot

perGranCorr :: CommonBlock -> MatrixCmplx -> MatrixReal -> RlxRoot -> MatrixCmplx
perGranCorr cb aMat vMat rlxRt = 
  let
     states       = getStates cb
     eTot         = getETot cb
     deco         = getDeco cb
     nSubStep     = getSubStep cb
     dt           = getDt cb
     rightIndex (n,m) = n*states + m
     rightTupla x = (div x states, mod x states)
     rlxRootI     = rlxRt -1
     vRlx         = vMat !! (rightIndex(rlxRootI,rlxRootI))
     eKin         = eTot - vRlx
     tau          = fmap (\x -> abs(1/(vRlx-(vMat!! rightIndex(x,x))))*(1+deco/eKin)) [0..states-1]
     firstSel     = fmap rightIndex [(x,y)| x<-[0..states-1], y <-[0..states-1], x/=rlxRootI, y/=rlxRootI]
     rightDT      = dt / (fromIntegral nSubStep)
     newE (x,y)   = (aMat!!(rightIndex (x,y)))*((exp(-rightDT/tau!!x)):+0.0)*((exp(-rightDT/tau!!y)):+0.0)
     firstCorr    = fmap (\(i,x) -> if elem i firstSel then (newE $ rightTupla i) else x) $ zip [0..] aMat
     secondSel    = fmap rightIndex [(x,y)| x<-[0..states-1], y<-[0..states-1], x==rlxRootI, y==rlxRootI]
     otherStatI   = fmap rightIndex [(x,y)| x<-[0..states-1], y<-[0..states-1], x==y,x/=rlxRootI]
     populOsSum   = sum $ fmap (realPart) (fmap (firstCorr !!) otherStatI)
     popRlxRt     = (1.0 - populOsSum) :+ 0.0
     secondCorr   = fmap (\(i,x) -> if elem i secondSel then popRlxRt else x) $ zip [0..] firstCorr
     thirdSel     = fmap rightIndex [(x,y)| x<-[0..states-1], y<-[0..states-1], x==rlxRootI ||  y==rlxRootI, (x,y)/=(rlxRootI,rlxRootI)]
     sqrtApAc     = sqrt ((popRlxRt)/(aMat!!(rightIndex (rlxRootI,rlxRootI))))
     newE' (x,y)  = if x/=rlxRootI then (secondCorr!!((rightIndex (x,y))))*((exp(-rightDT/tau!!x)):+0.0)*sqrtApAc 
                                   else (secondCorr!!((rightIndex (x,y))))*((exp(-rightDT/tau!!y)):+0.0)*sqrtApAc
     thirdCorr    = fmap (\(i,x) -> if elem i thirdSel then (newE' $ rightTupla i) else x) $ zip [0..] secondCorr
  in thirdCorr

aleGenRan :: SR.RandomGen g => g -> [Double]
aleGenRan seed = SR.randomRs (0.0,1.0) seed

printWell' :: [[Double]] -> IO()
printWell' xs = putStrLn $ unlines (fmap unwords (fmap (fmap $ printf "%12.10E") xs))

printWell :: MatrixReal -> Int -> String
printWell xs states = unlines (fmap unwords (fmap (fmap $ printBetter . printf "%+12.10e") $ chunksOf states xs))

printBetter fStr = let
   existChars char fStr = elem char (drop ((maybe (-1) id $ elemIndex 'e' fStr)+1) fStr)
   correct fStr         = (take ((maybe (-1) id $ elemIndex 'e' fStr)+1) fStr) ++ "+" ++ (drop ((maybe (-1) id $ elemIndex 'e' fStr)+1) fStr)
   tamp                 = if (existChars '+' fStr || existChars '-' fStr) then fStr else correct fStr
   in if (head tamp == '+') then " " ++ tail tamp else tamp

printComplex :: (RealFloat a, PrintfArg a) => Complex a -> String
printComplex fStr = let
           segno = if imagPart fStr > 0 then " +" else " -"
           in "(" ++ (printBetter $ printf "%+8.5e" (realPart fStr)) ++ segno ++ (printBetter $ printf "%+8.5e" (abs $ imagPart fStr)) ++ "i)"
   
printWellC :: MatrixCmplx -> Int -> String
printWellC xs states = unlines $ fmap unwords $ fmap (fmap printComplex) $ chunksOf states xs


--handshake :: (Eq a) => [a] -> [a] -> [[[a]]]
--handshake l s = let n = length l
--                    cartesian =[[x,y] | x <- l, y <- s]
--                    indexPermutation = permutations $ take n [0..]
--                    indexBase = repeat $ take n $iterate (+n) 0
--                in fmap (fmap (cartesian !!)) $ L.zipWith (L.zipWith(+)) indexPermutation indexBase
--
--scaleP :: Num a => [[a]] -> a               
--scaleP g = sum $ L.zipWith (*) (g!!0)  (g!!1)
--
--sumScalProd :: (Eq a, Num a) => [[a]] -> [[a]] -> [a]
--sumScalProd l s = fmap sum $ fmap (fmap scaleP) (handshake l s)
--
--sumScalProd2 :: (Eq a, Num a) => [[a]] -> [[a]] -> [a]
--sumScalProd2 l s = fmap (sum . fmap scaleP) $ handshake l s
--
--goodCouples :: (Num a, Ord a) => [[a]] -> [[a]] -> [[[a]]]
--goodCouples l s =  handshake l s !! (maybe (-1) id (elemIndex (maximum $ sumScalProd l s) (sumScalProd l s)))
--
--correctSigns :: (Num a, Ord a) => [[a]] -> [[a]] -> [[a]]
--correctSigns l s = L.zipWith (\item list -> fmap (item*) list) (fmap (signum . scaleP) (goodCouples l s)) s 

