import FFTW, LinearAlgebra, Integrals, DifferentialEquations
using RadiiPolynomial, GLMakie, TickTock

include("helpers.jl")
include("float_functions.jl")
include("approx_derivatives.jl")

__save__ = false
__plot__ = true
__save_location__ = "./figures/"

# Number of fourier coefficients we want
N::Int64 = 50
# Size of pertubation
ε::Float64 = 0
# Computing projection into higher space for FFT
# N_fft::Int64 = nextpow(2, 2*N+1)
N_fft::Int64 = 2^14

# Defining the space we are working in
ℱ = Fourier(N,1.0)
ℱ_pad = Fourier(div(N_fft,2),1.0)

# Building ellipse solution to Kepler problem

# Angles of rotation
# ψ = [rotation about x axis, rotation about y-axis, rotation about z-axis]
ψ::Vector{Float64} = [1;2;3]
e::Float64 = 1/2
ϕ::Float64 = 0

sample_points, sample_time = kepler_sample(ψ, e, ϕ, N_fft)

u = zeros(ComplexF64, ℱ^3)
component(u,1).coefficients[:] = 1/N_fft * project(Sequence(ℱ_pad, [FFTW.fftshift(FFTW.fft(sample_points[1,:])); 0]), ℱ).coefficients[:]
component(u,2).coefficients[:] = 1/N_fft * project(Sequence(ℱ_pad, [FFTW.fftshift(FFTW.fft(sample_points[2,:])); 0]), ℱ).coefficients[:]
component(u,3).coefficients[:] = 1/N_fft * project(Sequence(ℱ_pad, [FFTW.fftshift(FFTW.fft(sample_points[3,:])); 0]), ℱ).coefficients[:]

F = zeros(ComplexF64, ℱ^3)
DF = zeros(ComplexF64, ℱ^3, ℱ^3)
DF_approx = zeros(ComplexF64, ℱ^3, ℱ^3)

println("Evaluating function and its derivative")

# tick()
F!(F, u, ε, N_fft)
# tock()

# tick()
DF!(DF, u, ε, N_fft)
# tock()

# tick()
# DF_approx!(DF_approx, u, ε, N_fft)
# tock()

println("Size of kernel before Newton = " * string(size(LinearAlgebra.nullspace(DF.coefficients),2)) * "\n")
# println("Difference between true derivative and finite differences: " * string(opnorm(DF-DF_approx)) * "\n")

# Applying Newton's method we solve for a numerical solution
Newton!(u, ε, F, DF)

println("|u₁[N]| = " * string(norm(component(u,1)[N])) * ", |u₁[-N]| = " * string(norm(component(u,1)[-N])))
println("|u₂[N]| = " * string(norm(component(u,2)[N])) * ", |u₂[-N]| = " * string(norm(component(u,2)[-N])))
println("|u₃[N]| = " * string(norm(component(u,3)[N])) * ", |u₃[-N]| = " * string(norm(component(u,3)[-N])))

println("Size of kernel after Newton = " * string(size(LinearAlgebra.nullspace(DF.coefficients),2)) * "\n")

# Plotting using GLMakie
time_data = collect(LinRange(0, 2*pi, 1000))
u_data = collection_eval(time_data, u)

sol_plot = Figure(size = (1000, 600))
sol_ax = Axis3(sol_plot[1,1], title = L"\text{Approximate solution to the Kepler problem}", 
    titlesize = 20,
    xlabel = L"$u_1$",
    xlabelsize = 20,
    ylabel = L"$u_2$",
    ylabelsize = 20,
    zlabel = L"$u_3$",
    zlabelsize = 20,
    limits = ((-2, 2), (-2, 2), (-2, 2)) 
    )

GLMakie.lines!(sol_ax,
    u_data[1,:],
    u_data[2,:],
    u_data[3,:],
    label = L"\text{Periodic Solution}"
    )

GLMakie.scatter!(sol_ax,
    0,
    0,
    0,
    label = L"\text{Origin}"
    )

GLMakie.scatter!(sol_ax,
    sample_points[1,:],
    sample_points[2,:],
    sample_points[3,:],
    label = L"\text{Sample points}"
    )

axislegend("Legend")

coordinates_plot = Figure(size = (1600, 600))
x_ax = Axis(coordinates_plot[1,1], title = L"\text{$x$-coordinate of approximate solution}", 
    titlesize = 20,
    xlabel = L"$t$",
    xlabelsize = 20,
    ylabel = L"$x$",
    ylabelsize = 20,
    limits = ((0, 2*pi), (-2, 2)) 
    )

GLMakie.lines!(x_ax,
    time_data,
    u_data[1,:],
    label = L"\text{Periodic solution}"
    )

GLMakie.lines!(x_ax,
    time_data,
    zeros(1000,1)[:],
    label = L"\text{Origin}"
    )

GLMakie.scatter!(x_ax,
    sample_time,
    sample_points[1,:],
    label = L"\text{Sample points}"
    )

axislegend("Legend")

y_ax = Axis(coordinates_plot[1,2], title = L"\text{$y$-coordinate of approximate solution}", 
    titlesize = 20,
    xlabel = L"$t$",
    xlabelsize = 20,
    ylabel = L"$y$",
    ylabelsize = 20,
    limits = ((0, 2*pi), (-2, 2)) 
    )

GLMakie.lines!(y_ax,
    time_data,
    u_data[2,:],
    label = L"\text{Periodic solution}"
    )

GLMakie.lines!(y_ax,
    time_data,
    zeros(1000,1)[:],
    label = L"\text{Origin}"
    )

GLMakie.scatter!(y_ax,
    sample_time,
    sample_points[2,:],
    label = L"\text{Sample points}"
    )

axislegend("Legend")

z_ax = Axis(coordinates_plot[1,3], title = L"\text{$z$-coordinate of approximate solution}", 
    titlesize = 20,
    xlabel = L"$t$",
    xlabelsize = 20,
    ylabel = L"$z$",
    ylabelsize = 20,
    limits = ((0, 2*pi), (-2, 2)) 
    )

GLMakie.lines!(z_ax,
    time_data,
    u_data[3,:],
    label = L"\text{Periodic solution}"
    )

GLMakie.lines!(z_ax,
    time_data,
    zeros(1000,1)[:],
    label = L"\text{Origin}"
    )

GLMakie.scatter!(z_ax,
    sample_time,
    sample_points[3,:],
    label = L"\text{Sample points}"
    )

axislegend("Legend")

if __plot__
    display(GLMakie.Screen(), sol_plot)
    display(GLMakie.Screen(), coordinates_plot)
end

if __save__
    save(string(__save_location__, "ellipse_sol" * string(e) * ".png"), sol_plot, px_per_unit = 8)
    save(string(__save_location__, "coordinates_sol" * string(e) * ".png"), coordinates_plot, px_per_unit = 8)
end