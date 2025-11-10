import FFTW
using RadiiPolynomial, GLMakie

include("helpers.jl")
include("float_functions.jl")
include("approx_derivatives.jl")

# Number of fourier coefficients we want
N = 10
# Computing projection into higher space for FFT
N_fft = 2^10

# Defining the space we are working in
ℱ = Fourier(N,1.0)

S = ℱ^3

# Building u corresponding to the unit circle in R³
u = zeros(ComplexF64, S)
component(u,1)[1] = 1
component(u,1)[-1] = 1

component(u,2)[1] = 1im
component(u,2)[-1] = -1im

component(u,3)[2] = 1
component(u,3)[-2] = 1

#u = Sequence(S, rand(3*(2*N+1)))

F = zeros(ComplexF64, S)
DF = zeros(ComplexF64, S, S)
DF_approx = zeros(ComplexF64, S, S)

# DF!(DF, u, N_fft)
# DF_approx!(DF_approx, u, N_fft)

# println(maximum(abs.((DF-DF_approx).coefficients[:])))

# Applying Newton's method we solve for a numerical solution
u = Newton(u, F, DF,1e-12,100)

# Plotting using GLMakie
time_data = collect(LinRange(0, 2*pi, 1000))
u_data = collection_eval(time_data, u)

sol_plot = Figure()
sol_ax = Axis3(sol_plot[1,1], title = L"\text{Approximate solution to the Kepler problem}", 
    titlesize = 20,
    xlabel = L"$u_1$",
    xlabelsize = 20,
    ylabel = L"$u_2$",
    ylabelsize = 20,
    zlabel = L"$u_3$",
    zlabelsize = 20
    )

GLMakie.lines!(sol_ax,
    u_data[1,:],
    u_data[2,:],
    u_data[3,:],
    label = L"Periodic Solution"
    )

display(GLMakie.Screen(), sol_plot)








