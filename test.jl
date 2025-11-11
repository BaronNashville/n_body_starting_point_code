include("helpers.jl")

N = 2
S = Fourier(N, 1.0)

# N_fft = nextpow(2, 2*N+1)
N_fft = 2^4

S_pad = Fourier(div(N_fft,2), 1.0)

G(u) = u/norm(u)

# u = zeros(S^3)
# component(u,1)[:] = [-2;-1;0;1;2]
# component(u,2)[:] = [-4;-1;2;1;5]
# component(u,3)[:] = [3;-1;6;1;4]

u = Sequence(S, [-2;-1;2;1;2])

G_of_u = zeros(ComplexF64, S)
G_of_u!(G_of_u, u, G, N_fft)


