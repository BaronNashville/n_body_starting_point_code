#############################################################
# This script rigorously solves for θ in the solution       #
# to the Kepler problem θ = 1 / c^3 (1 + e cosθ)^2          #
#############################################################
include("fourier_computations.jl")
include("kepler_maps.jl")

# Fix eccentricity
e::Float64 = 0.5

# Setup the space on which (θ, c) lives
Ncheb = 100
space = Chebyshev(Ncheb) × ParameterSpace()

X = zeros(space)
