[![nimble](https://raw.githubusercontent.com/yglukhov/nimble-tag/master/nimble.png)](https://github.com/yglukhov/nimble-tag)
[![MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

# What is randgen?
A random variable generating library for nim. 

# Installation
```
git clone https://github.com/YesDrX/randgen.git
cd randgen
nimble install
```
or
```
nimble install randgen
```

# [Documentaion](https://yesdrx.github.io/randgen/)

# Example
```nim
import randgen
when isMainModule:
  discard initRand(100)
  echo randGaussian()
```

```
0.754597695851647
```
