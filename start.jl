# Setup the correct julia environment
import Pkg
Pkg.activate("nbody")
Pkg.resolve()

# Launch the program
#include("main.jl")

import FFTW, LinearAlgebra, Integrals, DifferentialEquations
using RadiiPolynomial, GLMakie, TickTock