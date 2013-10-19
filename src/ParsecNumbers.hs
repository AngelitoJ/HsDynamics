-- HsParser: A Parsec builder, a toy for experimenting things:
-- @2013 Angel Alvarez, Felipe Zapata, from The ResMol Group  


module ParsecNumbers where

import Data.Char
import Data.Functor.Identity
import Data.Either
import Control.Monad
import Text.Parsec
import Text.Parsec.Char
import Text.Parsec.Combinator

-- ========> Internal Modules
import CommonTypes

-- =================> TYPES <================
data Sign  = Neg | Pos deriving Show

-- Some silly numeric parser, (it may be a library in the hackage, whatsoever, we duplicated that here, just for fun)
-- Parse a real as <Integer>'.'<Decimal>['E'<Integer>] 
-- Usually : "['+','-']"[0-9]"."{8,[0-9]}"E"["-","+"][0-9][0-9]
realNumber :: MyParser st Double
realNumber = try $ do
    (integerPart,s)    <- integerNumberAndSign
    fractionalPart <- fractionNumber
    exponentPart   <- option 1.0 exponentNumber
    let baseNumber = case s of 
                       Neg -> negate $ (fromInteger integerPart) + fractionalPart
                       Pos -> (fromInteger integerPart) + fractionalPart
                                   
        result = case exponentPart of
                      1.0   -> baseNumber
                      other -> baseNumber * exponentPart
    return result

-- Parse a fractional part , including leading '.'
fractionNumber :: MyParser st Double
fractionNumber = do
    char '.'                                   -- leading dot
    digits <- many1 digit
    let l = length digits 
        n = fromIntegral $ foldl (\x d -> 10*x + toInteger (digitToInt d)) 0 digits
        result = seq n (n / 10^l)
    return result

-- Parse an exponent part 
exponentNumber :: MyParser st Double
exponentNumber = do
    oneOf "eE"
    f          <- sign
    e          <- decimalNumber
    let result = power (f e)
    seq result (return result)
    where
        power e  | e < 0      = 1.0/power(-e)
                 | otherwise  = fromInteger (10^e)


-- Parse an int as "-/+"{*,[0-9]} 
intNumber :: MyParser st Int
intNumber = do
    number <- integerNumber
    return $ fromIntegral number

-- Parse an integer as "-/+"{*,[0-9]} 
integerNumber :: MyParser st Integer
integerNumber = do
    f <- sign
    n <- decimalNumber
    return (f n)
    
integerNumberAndSign :: MyParser st (Integer,Sign)
integerNumberAndSign = do
    f <- sign'
    n <- decimalNumber
    return (n,f)
   
-- Parse sign '+-'
sign :: MyParser st (Integer -> Integer)
sign =
    (char '-' >> return negate) 
    <|> (char '+' >> return id)
    <|> return id
    
-- Parse sign '+-'
sign' :: MyParser st (Sign)
sign' =
    (char '-' >> return Neg) 
    <|> (char '+' >> return Pos)
    <|> return Pos

-- Parse a number in base 16 (computer geeks like this)
hexadecimalNumber :: MyParser st Integer
hexadecimalNumber = oneOf "xX" >> naturalNumber 16 hexDigit

-- Parse a number in base 10
decimalNumber :: MyParser st Integer
decimalNumber = naturalNumber 10 digit

-- Parse one or more digits into a suitable natural number (as an integer in base BaseDigit)
naturalNumber :: Integer -> MyParser st Char -> MyParser st Integer
naturalNumber base baseDigit = do
    digits <- many1 baseDigit
    let n = foldl (\x d -> base*x + toInteger (digitToInt d)) 0 digits
    seq n (return n)


-- No QuickCheck right now, so we resort to a poor's man test routine...
test = zip numbers $ rights $ map reals numbers
    where
        reals = parse realNumber ""
        numbers = [
             "-5.458429","-0.877245","0.459260"
            ,"-2.956134","-1.597542","-0.054604"
            ,"-1.041305","0.278925" ,"-0.248875"
            ,"1.253079" ,"-0.712285","-0.589929"
            ,"3.604078" ,"0.489114" ,"-0.043248"
            ,"5.668394" ,"-0.770555","0.360499"
            ,"1.154010" ,"-2.474177","-1.298762"
            ,"-1.731324","2.952315" ,"-0.080602"
            ,"3.561637" ,"2.563500" ,"0.619739"
            ,"-2.513848","-3.554577","0.058753"
            ,"7.271972" ,"0.058884" ,"1.018413"
            ,"5.772123" ,"-2.552650","0.638705"
            ,"-7.076648","-1.863584","0.723654"
            ,"-5.832434","1.004868" ,"0.733958"
            ,"0.048577" ,"4.073074" ,"-0.148143"
            ,"-2.656692","3.737933" ,"-2.051326"
            ,"-3.116744","3.363378" ,"1.340101"
            ]
