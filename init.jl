# Setup the correct julia environment
import Pkg
Pkg.activate("nbody")
Pkg.instantiate()
# Pkg.resolve()

import FFTW, LinearAlgebra, Integrals, DifferentialEquations
using RadiiPolynomial, LinearAlgebra, GLMakie, TickTock, DifferentialEquations

# Launch the program
#include("main.jl")
#include("test.jl")