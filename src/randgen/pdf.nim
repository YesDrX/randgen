import math

proc pdfGaussian*(mu = 0.0, sigma = 1.0, x = 0.0): float=
  ## Wiki: https://en.wikipedia.org/wiki/Normal_distribution
  ## 
  result = 1.0/(sigma*sqrt(2.0 * PI)) * exp(-0.5 * (x - mu) * (x-mu) / sigma)

when isMainModule:
  echo pdfGaussian()