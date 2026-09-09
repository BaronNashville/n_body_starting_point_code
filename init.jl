# Setup the correct julia environment
import Pkg
Pkg.activate("nbody")
Pkg.instantiate()
# Pkg.resolve()

import FFTW, LinearAlgebra, Integrals, DifferentialEquations, GeometryBasics
using RadiiPolynomial, LinearAlgebra, GLMakie, TickTock, DifferentialEquations, Colors, CSV, DataFrames

# Launch the program
#include("main.jl")
#include("test.jl")
# include("starting_point.jl")