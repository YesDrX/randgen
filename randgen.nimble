# Package

version       = "0.1.0"
author        = "YesDrX"
description   = "A random variable generator library for nim-lang."
license       = "MIT"
srcDir        = "src"



# Dependencies

requires "nim >= 1.0.0"

task donothing, "donothing":
    exec "nim c -r src/randgen.nim"
    rmFile "src/randgen"
