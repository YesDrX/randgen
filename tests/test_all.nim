import ../src/randgen
import random
import stats
import math
import strformat
import nimpy
import sequtils
import os

let
  py = pyBuiltinsModule()
  plt = pyImport("matplotlib.pyplot")
  sbn = pyImport("seaborn")

proc plot(data: seq[float], figname_short: string, figname_long: string)=
  discard plt.close()
  var
    fig_ax = plt.subplots().to(seq[PyObject])
    fig = fig_ax[0]
    ax = fig_ax[1]
  discard sbn.distplot(data, ax=ax)
  discard plt.title(fmt"Density for {figname_long}")
  # discard fig.savefig(absolutePath("tests/" & figname & ".png"))
  discard fig.savefig(absolutePath("./tests/" & figname_short & ".svg"))

proc  run()=
  randomize()
  var
    data_len = 100_000
    data_float = newSeq[float](data_len)
    diff_mean = 0.0
    diff_std = 0.0
    epsilon = 0.01
    expected_mean =  0.0
    expected_std = 0.0
    expected_skw = 0.0
    expected_kur = 0.0
    multiplier = 1.0

  # uniform
  for i in 0 ..< data_len:
    data_float[i] = rand(1.0)
  expected_mean = 0.5
  expected_std = 1/sqrt(12.0)
  expected_skw = 0.0
  expected_kur = -1.2
  multiplier = 1.5
  plot(data_float, "uniform", "Sample of Uniform[0,1]")

  doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier
  doAssert (data_float.kurtosis - expected_kur).abs < epsilon  * multiplier * multiplier * multiplier
  echo "Test 1:"
  echo fmt"Uniform[0,1]: sample mean/std/skewness/kurtosis passed! (epsilon =  {epsilon * multiplier})"

  # normal
  for i in 0 ..< data_len:
    data_float[i] = randGaussian(0.5, 0.9)
  expected_mean = 0.5
  expected_std = 0.9
  expected_skw = 0.0
  expected_kur = 0.0
  multiplier = 2.0
  plot(data_float, "normal", "Sample of Gaussian(0.5,0.9^2)")

  doAssert (data_float.mean() - expected_mean).abs < epsilon *  multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon *  multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon *  multiplier * multiplier * multiplier
  doAssert (data_float.kurtosis - expected_kur).abs < epsilon *  multiplier * multiplier * multiplier * multiplier
  echo "Test 2:"
  echo fmt"Normal[0.5,0.9]: sample mean/std/skewness/kurtosis passed! (epsilon =  {epsilon * multiplier})"

  # exponential
  for i in 0 ..< data_len:
    data_float[i] = randExponential(1.5)
  expected_mean = 1 / 1.5
  expected_std = 1 / 1.5
  expected_skw = 2.0
  expected_kur = 6.0
  multiplier = 3.0
  plot(data_float, "exponential", "Sample of Exponential(1.5)")

  doAssert (data_float.mean() - expected_mean).abs < epsilon *  multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon *  multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon *  multiplier * multiplier * multiplier
  doAssert (data_float.kurtosis - expected_kur).abs < epsilon * multiplier * multiplier * multiplier * multiplier
  echo "Test 3:"
  echo fmt"Exponential[1.5]: sample mean/std/skewness/kurtosis passed! (epsilon =  {epsilon * multiplier})"

  # laplace
  for i in 0 ..< data_len:
    data_float[i] = randLaplace(5.0, 4.0)
  expected_mean = 5.0
  expected_std = sqrt(2.0) * 4.0
  expected_skw = 0.0
  expected_kur = sqrt(3.0)
  multiplier = 5.0
  plot(data_float, "laplace", "Sample of Laplace(5.0,4.0)")

  doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier * multiplier
  doAssert (data_float.kurtosis - expected_kur).abs < epsilon * multiplier * multiplier * multiplier * multiplier * multiplier
  echo "Test 4:"
  echo fmt"Laplace[5,4]: sample mean/std/skewness/kurtosis passed! (epsilon =  {epsilon * multiplier})"

  # rayleigh
  for i in 0 ..< data_len:
    data_float[i] = randRayleigh(2.0)
  expected_mean = sqrt(PI / 2.0) * 2.0
  expected_std = sqrt((4.0 - PI) / 2.0) * 2.0
  expected_skw = 2.0 * sqrt(PI) * (PI - 3.0) / pow(4.0 - PI, 1.5)
  expected_kur = -( 6.0 * PI * PI - 24.0 * PI + 16.0) / pow(4.0 - PI, 2.0)
  multiplier = 5.0
  plot(data_float, "rayleigh", "Sample of Rayleigh(2.0)")

  doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier * multiplier
  doAssert (data_float.kurtosis - expected_kur).abs < epsilon * multiplier * multiplier * multiplier * multiplier * multiplier
  echo "Test 5:"
  echo fmt"Rayleigh[2.0]: sample mean/std/skewness/kurtosis passed! (epsilon =  {epsilon * multiplier})"

  # gamma
  for i in 0 ..< data_len:
    data_float[i] = randGamma(3.0, 2.0)
  expected_mean = 2.0 * 3.0
  expected_std = sqrt(3.0) * 2.0
  expected_skw = 2.0 / sqrt(3.0)
  expected_kur = 6.0 / 3.0
  multiplier = 3.0
  plot(data_float, "gamma", "Sample of Gamma(3.0,2.0)")

  doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier
  doAssert (data_float.kurtosis - expected_kur).abs < epsilon * multiplier * multiplier * multiplier * multiplier * multiplier
  echo "Test 6:"
  echo fmt"Gamma[3.0,2.0]: sample mean/std/skewness/kurtosis passed! (epsilon =  {epsilon * multiplier})"

  # lognormal
  for i in 0 ..< data_len:
    data_float[i] = randLogNormal(0.5, 0.9)
  plot(data_float, "lognormal", "Sample of Lognormal(0.5,0.9^2)")

  var
    mu = 0.5
    sigma = 0.9
  expected_mean = exp(mu + sigma * sigma / 2.0)
  expected_std = sqrt(exp(sigma * sigma - 1.0) * exp(2.0 * mu + sigma * sigma))
  expected_skw = (exp(sigma * sigma) + 2.0) * sqrt(exp(sigma * sigma) - 1.0)
  expected_kur = exp(4.0 * sigma * sigma) + 2.0 * exp(3.0 * sigma * sigma) + 3.0 * exp(2.0 * sigma * sigma) - 6.0
  multiplier = 100.0
  doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier * multiplier
  doAssert (data_float.kurtosis - expected_kur).abs < epsilon * multiplier * multiplier * multiplier * multiplier * multiplier
  echo "Test 7:"
  echo fmt"LogNormal[0.5,0.9]: sample mean/std/skewness/kurtosis passed! (epsilon =  {epsilon * multiplier})"

  # chi2
  for i in 0 ..< data_len:
    data_float[i] = randChi2(4)
  expected_mean = 4.0
  expected_std = sqrt(2.0  * 4.0)
  expected_skw = sqrt(8.0 / 4.0)
  expected_kur = 12.0 / 4.0
  multiplier = 2.5
  plot(data_float, "chi2", "Sample of Chi2(4)")

  doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier * multiplier
  doAssert (data_float.kurtosis - expected_kur).abs < epsilon * multiplier * multiplier * multiplier * multiplier * multiplier
  echo "Test 8:"
  echo fmt"Chi2[4]: sample mean/std/skewness/kurtosis passed! (epsilon =  {epsilon * multiplier})"

  # f
  for i in 0 ..< data_len:
    data_float[i] = randF(7,10)
  var
    d1 = 7.0
    d2 = 10.0
  plot(data_float, "f", "Sample of F-Distribution(7,10)")

  expected_mean = d2 / (d2 - 2.0)
  expected_std = sqrt((2.0*d2*d2*(d1+d2-2.0))/(d1*(d2-2.0)*(d2-2.0)*(d2 - 4.0)))
  expected_skw = (2*d1+d2-2.0)*sqrt(8.0*(d2-4.0))/((d2-6.0)*sqrt(d1*(d1+d2-2.0)))
  multiplier = 3
  doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier * multiplier
  echo "Test 9:"
  echo fmt"F[7,10]: sample mean/std/skewness passed! (epsilon =  {epsilon * multiplier})"

  # t
  for i in 0 ..< data_len:
    data_float[i] = randStudent(5)
  var
    df = 5.0
  expected_mean = 0.0
  expected_std = sqrt(df / (df - 2.0))
  expected_skw = 0.0
  expected_kur = 6.0 / (df - 4.0)
  multiplier = 10.0
  plot(data_float, "student", "Sample of Student's t (4)")
  doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier * multiplier
  echo "Test 9:"
  echo fmt"student[5]: sample mean/std/skewness passed! (epsilon =  {epsilon * multiplier})"

  # betas
  for i in 0 ..< data_len:
    data_float[i] = randBeta(2.0,5.0)
  plot(data_float, "beta", "Sample of Beta(2.0,5.0)")
  var
    alpha = 2.0
    beta = 5.0
  expected_mean = alpha / (alpha + beta)
  expected_std = sqrt(alpha * beta / ((alpha + beta) * (alpha + beta) * (alpha + beta + 1.0)))
  expected_skw = 2 * (beta - alpha) * sqrt(alpha + beta + 1.0)/((alpha + beta + 2.0) * sqrt(alpha * beta))
  expected_kur = 6 * ((alpha - beta)^2 * (alpha + beta + 1.0) - alpha * beta * (alpha + beta + 2.0)) / (alpha * beta * (alpha + beta + 2.0)*(alpha + beta + 3.0))
  multiplier = 3.0
  doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier * multiplier
  doAssert (data_float.kurtosis - expected_kur).abs < epsilon * multiplier * multiplier * multiplier * multiplier * multiplier
  echo "Test 10:"
  echo fmt"beta[2.0,5.0]: sample mean/std/skewness passed! (epsilon =  {epsilon * multiplier})"

  # logistic
  for i in 0 ..< data_len:
    data_float[i] = randLogistic(9.0, 4.0)
  plot(data_float, "logistic", "Sample of Logistic(9.0,4.0)")
  var
    s = 4.0
  mu = 9.0
  expected_mean = mu
  expected_std = s * PI / sqrt(3.0)
  expected_skw = 0.0
  expected_kur = 1.2
  multiplier = 5.0
  doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier * multiplier
  doAssert (data_float.kurtosis - expected_kur).abs < epsilon * multiplier * multiplier * multiplier * multiplier * multiplier
  echo "Test 11:"
  echo fmt"logistic[9.0,4.0]: sample mean/std/skewness/kurtosis passed! (epsilon =  {epsilon * multiplier})"

  # # weibull
  # for i in 0 ..< data_len:
  #   data_float[i] = randWeibull(1.0, 1.5)
  # plot(data_float, "weibull")
  # expected_mean = mu
  # expected_std = s * PI / sqrt(3.0)
  # expected_skw = 0.0
  # expected_kur = 1.2
  # multiplier = 3.0
  # doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  # doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  # doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier * multiplier
  # doAssert (data_float.kurtosis - expected_kur).abs < epsilon * multiplier * multiplier * multiplier * multiplier * multiplier
  # echo "Test 12:"
  # echo fmt"Weibull[9.0,4.0]: sample mean/std/skewness/kurtosis passed! (epsilon =  {epsilon * multiplier})"

  # binomial
  for i in 0 ..< data_len:
    data_float[i] = randBinomial(0.6, 100).float
  plot(data_float, "binomial", "Sample of Binomial(0.6, 100)")
  var
    p = 0.6
    q = 1-p
    n = 100.0
  expected_mean = n * p
  expected_std = sqrt(n * p * q)
  expected_skw = (q-p) / sqrt(n*p*q)
  expected_kur = (1 - 6.0 * p * q)/(n*p*q)
  multiplier = 60.0
  doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier * multiplier
  doAssert (data_float.kurtosis - expected_kur).abs < epsilon * multiplier * multiplier * multiplier * multiplier * multiplier
  echo "Test 13:"
  echo fmt"binomial[0.6,100]: sample mean/std/skewness/kurtosis passed! (epsilon =  {epsilon * multiplier})"

  # poisson
  for i in 0 ..< data_len:
    data_float[i] = randPoisson(5.0).float
  plot(data_float, "poisson", "Sample of Poisson(5.0)")
  var
    lam = 5.0
  expected_mean = lam
  expected_std = sqrt(lam)
  expected_skw = 1.0 / sqrt(lam)
  expected_kur = 1.0 / lam
  multiplier = 3.0
  doAssert (data_float.mean() - expected_mean).abs < epsilon * multiplier
  doAssert (data_float.standardDeviation - expected_std).abs < epsilon * multiplier
  doAssert (data_float.skewness - expected_skw).abs < epsilon * multiplier * multiplier * multiplier * multiplier
  doAssert (data_float.kurtosis - expected_kur).abs < epsilon * multiplier * multiplier * multiplier * multiplier * multiplier
  echo "Test 14:"
  echo fmt"poisson[5.0]: sample mean/std/skewness/kurtosis passed! (epsilon =  {epsilon * multiplier})"

run()