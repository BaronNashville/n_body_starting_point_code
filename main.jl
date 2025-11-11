import FFTW, LinearAlgebra, Integrals, DifferentialEquations
using RadiiPolynomial, GLMakie

include("helpers.jl")
include("float_functions.jl")
include("approx_derivatives.jl")

# Number of fourier coefficients we want
N = 10
# Computing projection into higher space for FFT
N_fft = nextpow(2, 2*N+1)

# Defining the space we are working in
ℱ = Fourier(N,1.0)
ℱ_pad = Fourier(div(N_fft,2),1.0)

# Building ellipse solution to Kepler problem

# Angles of rotation
# ψ = [rotation about x axis, rotation about y-axis, rotation about z-axis]
Ψ = [0.0;0.0;0.0]
e = 0.1
ϕ = 0.0

sample_points, sample_time = kepler_sample(Ψ, e, ϕ, N_fft)

u = zeros(ComplexF64, ℱ^3)
component(u,1).coefficients[:] = 1/N_fft * project(Sequence(ℱ_pad, [FFTW.fftshift(FFTW.fft(FFTW.ifftshift(sample_points[1,:]))); 0]), ℱ).coefficients[:]
component(u,2).coefficients[:] = 1/N_fft * project(Sequence(ℱ_pad, [FFTW.fftshift(FFTW.fft(FFTW.ifftshift(sample_points[2,:]))); 0]), ℱ).coefficients[:]
component(u,3).coefficients[:] = 1/N_fft * project(Sequence(ℱ_pad, [FFTW.fftshift(FFTW.fft(FFTW.ifftshift(sample_points[3,:]))); 0]), ℱ).coefficients[:]

println("|u₁[N]| = " * string(norm(component(u,1)[N])) * ", |u₁[-N]| = " * string(norm(component(u,1)[-N])))
println("|u₂[N]| = " * string(norm(component(u,2)[N])) * ", |u₂[-N]| = " * string(norm(component(u,2)[-N])))
println("|u₃[N]| = " * string(norm(component(u,3)[N])) * ", |u₃[-N]| = " * string(norm(component(u,3)[-N])))

F = zeros(ComplexF64, ℱ^3)
DF = zeros(ComplexF64, ℱ^3, ℱ^3)
DF_approx = zeros(ComplexF64, ℱ^3, ℱ^3)

F!(F, u, N_fft)
println(norm(F))

DF!(DF, u, N_fft)
display(abs.(LinearAlgebra.eigvals(DF.coefficients)))

# DF!(DF, u, N_fft)
# DF_approx!(DF_approx, u, N_fft)

# println(opnorm(DF-DF_approx))

# # Applying Newton's method we solve for a numerical solution
# u = Newton(u, F, DF)

# println(LinearAlgebra.eigvals(DF.coefficients))

# Plotting using GLMakie
time_data = collect(LinRange(-pi, pi, 1000))
u_data = collection_eval(time_data, u)

sol_plot = Figure()
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
display(GLMakie.Screen(), sol_plot)

coordinates_plot = Figure()
x_ax = Axis(coordinates_plot[1,1], title = L"\text{$x$-coordinate of approximate solution}", 
    titlesize = 20,
    xlabel = L"$t$",
    xlabelsize = 20,
    ylabel = L"$x$",
    ylabelsize = 20,
    limits = ((-pi, pi), (-2, 2)) 
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
    limits = ((-pi, pi), (-2, 2)) 
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
    limits = ((-pi, pi), (-2, 2)) 
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

display(GLMakie.Screen(), coordinates_plot)