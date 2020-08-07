import math

proc cdfGaussian*(mu = 0.0, sigma = 1.0, x = 0.0): float=
  ## Winitzki (2008)
  ## 
  ## A Handy Approximation for the Error Function and its Inverse (https://www.academia.edu/9730974/A_handy_approximation_for_the_error_function_and_its_inverse)
  ## 
  var
    z = (x - mu) / sigma
  result = 0.5 * (1 + sqrt(1 - exp((-0.5*z*z*(4.0/PI + 0.147*0.5*z*z )/(1+0.147*0.5*z*z)))))

when isMainModule:
  echo cdfGaussian(x = 1.6)