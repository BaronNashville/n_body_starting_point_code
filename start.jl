# Setup the correct julia environment
import Pkg
Pkg.activate("nbody")
Pkg.instantiate()

import FFTW
using RadiiPolynomial
# Launch the program
# include("main.jl")