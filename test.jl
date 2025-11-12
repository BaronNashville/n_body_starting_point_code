x = 1

println(x+5)

# include("helpers.jl")
# include("float_functions.jl")
# include("approx_derivatives.jl")

# N = 1
# S = Fourier(N, 1.0)

# N_fft = nextpow(2, 2*N+1)

# S_pad = Fourier(div(N_fft,2), 1.0)

# G(u) = u/norm(u)^3

# # u = Sequence(S^3, rand(3*(2*N+1)))
# u = Sequence(S^3, [1;2;3;4;5;6;7;8;9])

# F = zeros(ComplexF64, S^3)
# DF = zeros(ComplexF64, S^3, S^3)
# DF_approx = zeros(ComplexF64, S^3, S^3)

# F!(F, u, N_fft)
# #DF!(DF, u, N_fft)
# #DF_approx!(DF_approx, u, N_fft)




