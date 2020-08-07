import random
import math

type
  FloatPair* = tuple[num1: float, num2: float]

proc randGaussian*(mu = 0.0, sigma = 1.0, return_pari : bool): FloatPair =
  ## Use Box-Muller method to generate a pair of Gaussian numbers.
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Box-Muller_transform
  ## 
  
  var
    u0 = rand(1.0)
    u1 = rand(1.0)
    z0 = sqrt( -2.0 * ln(u0) ) * cos( 2.0 * PI * u1)
    z1 = sqrt( -2.0 * ln(u1) ) * sin( 2.0 * PI * u0)
  z0 = z0 * sigma + mu
  z1 = z1 * sigma + mu
  result = (z0,z1)

proc randGaussian*(mu = 0.0, sigma = 1.0): float =
  ## Use Box-Muller method to generate a pair of Gaussian numbers.
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Box-Muller_transform
  ## 
  ## <img src="normal.svg" alt="Italian Trulli">
  ## 
    
  var
    (z0, z1) = randGaussian(mu, sigma, true)
  result = z0

proc randGaussianTrunc*(mu = 0.0, sigma = 1.0, lb = -2.0, ub = 2.0): float =
  var
    u = int.low.float  
  while not (u >= lb and u <= ub):
    u = randGaussian(mu, sigma)
  return u

proc randGaussianTanh*(mu = 0.0, sigma = 1.0, lb = -2.0, ub = 2.0): float =
  var
    u = randGaussian(0.0, 1.0).tanh
  return (u + 1.0) * 0.5 * (ub - lb) + lb

proc randExponential*(lambda = 1.0): float =
  ## Use inverse CDF method to generate a Exponetial random numer.
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Exponential_distribution
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Inverse_transform_sampling
  ## 
  ##  <img src="exponential.svg" alt="Italian Trulli">
  ## 

  var
    u = rand(1.0)
  result = -(ln(1.0-u)) / lambda

proc randLaplace*(mu = 0.0, b = 1.0): float =
  ## Use inverse CDF method to generate a Laplace random numer.
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Laplace_distribution
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Inverse_transform_sampling
  ## 
  ##  <img src="laplace.svg" alt="Italian Trulli">
  ## 
  
  var
    u = rand(1.0)
  if u >= 0.5:
    result = mu - b * ln( 2.0 * ( 1.0 - u))
  else:
    result = mu + b * ln( 2.0 * u)

proc randCauchy*(x0 = 0.0, gamma = 1.0): float =
  ## Use inverse CDF method to generate a Cauchy random numer.
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Cauchy_distribution
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Inverse_transform_sampling
  ## 

  var
    u = rand(1.0)
  result = x0 + gamma * tan( PI * (u - 0.5))

proc randRayleigh*(sigma = 1.0): float =
  ## Use inverse CDF method to generate a Rayleigh random numer.
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Rayleigh_distribution
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Inverse_transform_sampling
  ## 
  ##  <img src="rayleigh.svg" alt="Italian Trulli">
  ## 

  var
    u = rand(1.0)
  result = sqrt( -2.0 * ln(u) ) * sigma

proc randGammaLargeK(k = 2.0): float =
  ## Paper : Computer methods for sampling from gamma, beta, poisson and bionomial distributions
  ## 
  ## Ref : J. H. Ahrens and U. Dieter, Computing volume 12, pages223–246(1974)
  ## 
  ## Link : https://link.springer.com/article/10.1007/BF02293108
  ## 
  var
    sqa, x, y, v : float
  
  sqa = sqrt( 2.0 * k - 1.0)
  y = tan(PI * rand(1.0))
  x = sqa * y + k - 1.0
  v = rand(1.0)

  while (v > (1.0+y*y) * exp((k-1.0) * ln(x / (k - 1.0)) - sqa * y)):
    while (x < 0.0):
      y = tan(PI * rand(1.0))
      x = sqa * y + k - 1.0
    v = rand(1.0)
  result = x

proc randGammaFracK(k = 0.5): float=
  var
    p, q, x, u, v: float
  
  if k == 0.0:
    result = 0.0
  
  p = E / (k + E)
  u = rand(1.0)
  v = rand(1.0)
  if u < p:
    x = exp(1.0 / k * ln(v))
    q = exp(-x)
  else:
    x = 1 - ln(v)
    q = exp((k-1) * ln(x))
  
  while rand(1.0) >= q:
    p = E / (k + E)
    u = rand(1.0)
    v = rand(1.0)
    if u < p:
      x = exp(1.0 / k * ln(v))
      q = exp(-x)
    else:
      x = 1 - ln(v)
      q = exp((k-1) * ln(x))
  
  result = x

proc randGammaIntK(k = 5.0): float=
  if k < 12.0:
    var
      prod = 1.0
    for i in 0 ..< k.int:
      prod *= rand(1.0)
    result = -ln(prod)
  else:
    result = randGammaLargeK(k)

proc randGamma*(k = 1.0, theta = 2.0): float =
  ## Wiki: https://en.wikipedia.org/wiki/Gamma_distribution
  ## 
  ## <img src="gamma.svg" alt="Italian Trulli">
  ## 
  if k >= int32.high.float:
    result = theta * (randGammaLargeK(k.floor) + randGammaFracK(k - k.floor))
  elif k.floor == k:
    result = theta * (randGammaIntK(k.floor))
  elif k.floor == 0.0:
    result = theta * (randGammaFracK(k - k.floor))
  else:
    result = theta * (randGammaIntK(k.floor) + randGammaFracK(k - k.floor))

proc randLogNormal*(mu = 0.0, sigma = 1.0): float=
  ## Wiki: https://en.wikipedia.org/wiki/Log-normal_distribution
  ## 
  ## Method: first generate a Gaussian random number, and then take it to exp(.) function.
  ## 
  ## <img src="lognormal.svg" alt="Italian Trulli">
  ## 
  var
    rand_normal = randGaussian(mu, sigma)
  result = exp(rand_normal)

proc randChi2*(df = 1): float=
  ## Wiki: https://en.wikipedia.org/wiki/Chi-square_distribution
  ## 
  ## Method: The chi-square distribution (with d.o.f. df) is a special case of gamma distribution (Gamma(df/2, 2))
  ## 
  ## <img src="chi2.svg" alt="Italian Trulli">
  ## 
  result = 2.0 * randGamma(df.float / 2.0, 1.0)

proc randF*(df1 = 1, df2 = 1): float=
  ## Wiki: https://en.wikipedia.org/wiki/F-distribution
  ## 
  ## Method: F-distribution is defined as ratio of two independent and scaled chi-squared random variables.
  ## 
  ## <img src="f.svg" alt="Italian Trulli">
  ## 
  
  var
    y1 = randChi2(df1) / df1.float
    y2 = randChi2(df2) / df2.float
  result = y1 / y2

proc randStudent*(df = 1): float=
  ## Wiki: https://en.wikipedia.org/wiki/Student%27s_t-distribution
  ## 
  ## Method: Student t distribution is defined as ratio of standard normal random number devided by a Chi2 random numer (with d.o.f of df)
  ## 
  ## <img src="student.svg" alt="Italian Trulli">
  ## 
  if df <= 10:
    var
      y1 = randGaussian()
      y2 = sqrt(randChi2(df) / df.float)
    result = y1 / y2
  else:
    var
      y1, y2, z: float
    y1 = randGaussian()
    y2 = randExponential(1.0 / (df.float / 2.0 - 1.0))
    z = y1 * y1 / (df.float - 2.0)

    while 1.0 - z < 0.0 or 1.0 < exp(-y2 - z) + z:
      y1 = randGaussian()
      y2 = randExponential(1.0 / (df.float / 2.0 - 1.0))
      z = y1 * y1 / (df.float - 2.0)
    
    result = y1 / sqrt( (1.0 - 2.0 / df.float) * (1.0 - z))

proc randBeta*(alpha = 0.5, beta = 0.5): float=
  ## Wiki: https://en.wikipedia.org/wiki/Beta_distribution
  ## 
  ## <img src="beta.svg" alt="Italian Trulli">
  ##


  if alpha <= 1.0 and beta <= 1.0:
    var
      u,v,x,y : float
    while true:
      u = rand(1.0)
      v = rand(1.0)
      x = pow(u, 1.0 / alpha)
      y = pow(v, 1.0 / beta)
      if x + y <= 1.0:
        if x + y > 0.0:
          return x/(x+y)
        else:
          var
            logx = ln(u) / alpha
            logy = ln(v) / beta
            logm = if logx > logy : logx else: logy
          logx -= logm
          logy -= logm
          return exp(logx - ln(exp(logx) + exp(logy)))
  else:
    var
      x1 = randGamma(alpha, 1.0)
      x2 = randGamma(beta, 1.0)
    return x1 / (x1 + x2)

proc randLogistic*(mu = 5.0, s = 2.0): float=
  ## Use inverse CDF method to generate a Logistic random numer.
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Logistic_distribution
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Inverse_transform_sampling
  ## 
  ## <img src="logistic.svg" alt="Italian Trulli">
  ## 
  
  var
    u = rand(1.0)
  while u == 0.0:
    u = rand(1.0)
  result = mu - s * ln(1.0/u - 1.0)

proc randWeibull*(lambda = 1.0, k = 0.5): float=
  ## Use inverse CDF method to generate a Weibull random numer.
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Weibull_distribution
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Inverse_transform_sampling
  ## 
  var
    u = rand(1.0)
  while u == 1.0:
    u = rand(1.0)
  result = lambda * pow(-ln(1.0 - u), 1.0 / k)
  
proc randBinomial*(p = 0.5, n = 1): int=
  ## Use a truncated normal distribution to approximate the Binomial distribution.
  ## 
  ## Wiki: https://en.wikipedia.org/wiki/Binomial_distribution
  ## 
  ## <img src="binomial.svg" alt="Italian Trulli">
  ## 

  var
    q = 1-p
    k = 0
    u: float
  
  if n > 20:
    k = randGaussianTrunc((n-20).float*p,sqrt((n-20).float*p*q),0.0,(n-20).float).int
    # k = randGaussianTanh((n-20).float*p,sqrt((n-20).float*p*q),0.0,(n-20).float).int
  for i in 0 ..< min(20,n):
    u = rand(1.0)
    if u <= p: k += 1
  return k

proc randPoisson*(lambda = 1.0): int=
  ## Wiki: https://en.wikipedia.org/wiki/Poisson_distribution
  ## 
  ## <img src="poisson.svg" alt="Italian Trulli">
  ## 

  var
    elambda: float
    prod = 1.0
    m, x: float
    k = 0
    lam = lambda
  
  while lam > 10:
    m = lam * (7.0 / 8.0)
    x = randGamma(m)
    if x >= lam:
      return k.int + randBinomial(lam / x, m.int - 1)
    else:
      k += m.int
      lam -= x
  
  elambda = exp(-lam)
  while prod > elambda:
    prod *= rand(1.0)
    k += 1
  return k - 1  

when isMainModule:
  # discard initRand(100)
  randomize()
  # echo randGaussian(0.0, 1.0, true)
  # echo randGaussianTanh()
  ## echo randExponential()
  ## echo randLaplace()
  ## echo randCauchy()
  ## echo randRayleigh()
  # echo randGammaLargeK(100)
  # echo randGammaFracK(0.5)
  # echo randGamma()
  # echo randLogNormal()
  # echo randF()
  # echo randStudent()
  # echo randBeta()
  # echo randLogistic()
  # echo randWeibull()
  # echo randPoisson(100.0)
  echo randBinomial(0.6, 100)
  # echo randGaussianTrunc(0.0, 1.0, -2.0, 2.0)