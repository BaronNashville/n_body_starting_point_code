include("helpers.jl")
include("float_functions.jl")
include("approx_derivatives.jl")

# N = 1
# S = Fourier(N, 1.0)

# # N_fft = nextpow(2, 2*N+1)
N_fft = 2^14

# S_pad = Fourier(div(N_fft,2), 1.0)

# u = Sequence(S^3, rand(3*(2*N+1)))
# #u = Sequence(S^3, [1;2;3;4;5;6;7;8;9])

# F = zeros(ComplexF64, S^3)
# DF = zeros(ComplexF64, S^3, S^3)
# DF_approx = zeros(ComplexF64, S^3, S^3)

# F!(F, u, N_fft)
# DF!(DF, u, N_fft)
# DF_approx!(DF_approx, u, N_fft)
# println(opnorm(DF-DF_approx))

crit_points = [
    0.24        1.57    4.7     -0.48   3.1992146871613962
    -0.0090     4.71    4.7     4.8     0.28830893673407454
    0.75        4.9     6.5     2.8     0.4778866133081419
    0           0       0       0       3.573680206402789e-7
]

ψ::Vector{Float64} = 2*pi*rand(3)
e::Float64 = 0.1

X = [ψ;e]

X , success = newton((X) -> (DA_1(X, N_fft), HA_1_approx(X, N_fft)), X, maxiter = 20, verbose = true)

display(success)
display(X)





