-- CAM-ES algorithm for optimization
-- Version : 1.0
-- Reference : The CMA Evolution Strategy: A Comparing Review
--             Nikolaus Hansen

import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Util
import Data.List


data Parameter = Parameter {lambda :: Int,
                            mu :: Double,
                            w :: Vector Double,
                            c_sig :: Double,
                            d_sig :: Double,
                            c_c :: Double,
                            mu_cov :: Double,
                            c_cov :: Double} deriving Show

-- setting default parameter
default_parameter :: Int -> Parameter
default_parameter n = Parameter {lambda = lambda',
                                 mu = mu',
                                 w = w',
                                 c_sig = c_sig',
                                 d_sig = d_sig',
                                 c_c = c_c',
                                 mu_cov = mu_cov',
                                 c_cov = c_cov'}
                      where
                        lambda' = 4 + floor ( 3 * log ( fromIntegral n ) )
                        mu' = fromIntegral $ div lambda' 2 
                        tmp' = fromList [1.0..mu']
                        tmp = log tmp'
                        low = mu' * log(succ mu') - sumElements tmp
                        w' =  scale (1.0/low) ( constant (log (succ mu')) (div lambda' 2) - tmp )
                        mu_eff = 1.0 / (sumElements $ w'^2)
                        c_sig' = (mu_eff + 2.0) / (fromIntegral n + mu_eff + 3.0)
                        d_sig' = 1.0 + c_sig' + 2.0 * max 0 (sqrt ((mu_eff-1.0)/(fromIntegral n + 1.0)) - 1 )
                        c_c' = 4.0 / (fromIntegral n + 4)
                        mu_cov' = mu_eff
                        c_cov' = part1 + (1.0-1.0/mu_cov')* min 1.0 ((2.0*mu_eff-1.0)/((fromIntegral n + 2.0)^2 + mu_eff))
                        part1 =  1.0 / mu_cov' * 2.0 / (fromIntegral n + sqrt 2.0)^2

-- generate random vector from normal distribution 
normal :: Matrix Double -> Matrix Double -> Double -> Vector Double -> [Vector Double] -> [Vector Double]
normal b d sigma xmean arz = normal' b d sigma xmean arz []
normal' :: Matrix Double -> Matrix Double -> Double -> Vector Double -> [Vector Double] -> [Vector Double] -> [Vector Double]
normal' b d sigma xmean [] r = r
normal' b d sigma xmean (x:xs) r = normal' b d sigma xmean xs (r++[y])
  where
    z = b <> d <> x
    y = xmean + scale sigma z 

-- triu : takes the upper triangle of a matrix
aux1 :: [Double] -> Int -> [Double]
aux1 c 0 = c
aux1 c i = (toList $ constant 0 i) ++ drop i c
aux2 :: [Vector Double] -> Int -> [Vector Double] -> [Vector Double]
aux2 [] i r = r
aux2 (c:cs) i r =
  aux2 cs (i+1) (r ++ [fromList $aux1 (toList c) i])
triu :: Matrix Double -> Int -> Matrix Double
triu c n = fromRows $ aux2 (toRows c) n []

-- the evolution procedure
evolution :: (Vector Double -> Double) -> Int -> Int -> Int -> Double -> Int -> Vector Double -> Vector Double -> Matrix Double -> Double -> Vector Double -> Vector Double -> Parameter -> IO (Vector Double)

evolution feval counteval eigeneval stopeval stopfitness n p_sig p_c c sigma xmean x param =
  if counteval >= stopeval then return x
  else
    do     
      arz <- randn n $ lambda param
      let 
        (d',b) = eigSH c
        d = diag $ sqrt d'
        arx = normal b d sigma xmean $ toColumns arz
        arfitness = map feval arx
        (_,ind) = unzip $ sort $ zip arfitness [0..]
        mju  = div (lambda param) 2
        arindex = take mju ind
        xmean' = (trans $ extractRows arindex $ trans $ fromColumns arx) <> (w param)
        zmean' = (trans $ extractRows arindex $ trans arz) <> (w param)
        c_s = c_sig param
        mueff = mu_cov param
        p_sig_1' = scale (1.0 - c_s) p_sig
        p_sig_2' = sqrt (c_s * (2.0 - c_s) * mueff)
        p_sig_3' = b <> zmean'
        p_sig' = p_sig_1' + scale p_sig_2' p_sig_3'

        chiN = (fromIntegral n)**0.5*(1.0-1.0/(4.0*(fromIntegral n)) + 1.0/(21.0+(fromIntegral n)**2.0))
        lda = fromIntegral $ lambda param
        scala = sqrt (1.0 - (1.0 - c_s)**(2.0*(fromIntegral counteval)/lda))
        hsig' = norm p_sig' / scala / chiN - 1.5 - 1.0/(fromIntegral n - 0.5)
        hsig = if hsig' < 0 then 1 else 0
        cc = c_c param
        p_c_1' = scale (1.0-cc) p_c
        p_c_2' = hsig * sqrt (cc * (2.0-cc) * mueff)
        p_c_3' = b <> d <> zmean'
        p_c' = p_c_1' + scale p_c_2' p_c_3'

        ccov = c_cov param
        mucov = mu_cov param
        c_1' = scale (1.0-ccov) c
        c_2' = scale (ccov*(1.0/mucov)) c_3'
        c_3' = outer p_c' p_c' + c_4'
        c_4' = scale (cc * (2.0-cc) * (1.0-hsig))c
        c_5' = ccov * (1.0-1.0/mucov)
       
        weights = w param
        c_6' = b <> d <> (trans $ extractRows arindex $ trans arz)
        c_7' = diag weights <>  trans c_6'
        c_8' = scale c_5' $ c_6' <> c_7'
        c'' = c_1' + c_2' + c_8'

        cs = c_sig param
        damps =  1.0+2.0 * (max 0 (sqrt ((mueff-1.0)/(1.0+fromIntegral n)) - 1)) + cs
        sigma' = sigma * exp ((cs / damps) * (norm p_sig' / chiN - 1.0))

        sigma'' = if head arfitness == arfitness !! (min (1+floor(lda/2.0)) (2+ceiling(lda/4.0)))
                  then
                    sigma' * exp (0.2+cs/damps)
                  else
                    sigma'

        eigeneval' = if fromIntegral (counteval-eigeneval) > lda/ccov/(fromIntegral n)/10.0 then
                       counteval + 1
                     else
                       eigeneval
        c' = triu c 0 + (trans $ triu c 1)
      
      if head arfitness <= stopfitness then return $ arx !! (head arindex)
        else
        evolution feval (counteval+1) eigeneval' stopeval stopfitness n p_sig' p_c' c' sigma'' xmean' (arx !! (head arindex)) param
  

-- main CAM-ES function
camES :: (Vector Double -> Double) -> Int -> IO (Vector Double)
camES feval n = do
  xmeans' <- rand n 1
  evolution feval counteval eigeneval stopeval stopfitness n p_sig p_c c sigma (head $ toColumns xmeans') p_c param
  where
    eigeneval = 0
    p_sig = constant 0 n
    p_c = constant 0 n
    param = default_parameter n
    counteval = 0
    stopeval = n*n*1000
    stopfitness = 1.0e-10
    c = ident n
    sigma = 0.5

-- objective function
f :: Vector Double -> Double
f x = (x1-1.0)^2 + (x2-4.5)^2 + 10.0
  where x1 = x @> 0
        x2 = x @> 1
        
main = do
  x <- camES f 2
  (putStr . vecdisp (dispf 2)) x 
