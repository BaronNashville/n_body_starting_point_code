function G_of_u!(G_of_u::Sequence, u::Sequence, G::Function, N_fft::Int64)
    # INPUTS:
    # 1) u  = {u_n}_{n=-N}^{N} sequence of fourier coefficients
    # 2) g : R → R analytic function

    # Set ouput to 0
    G_of_u .= 0

    # Extract parameters from inputs
    S = space(u)
    f = frequency(S)

    S_pad = Fourier(div(N_fft,2), f)    

    u_pad = project(u, S_pad)

    # Evaluate the fourier series at grid points using the IFFT
    grid_points = Sequence(S_pad, [FFTW.fft(FFTW.ifftshift(u_pad.coefficients[begin:end-1])); 0])
    # Evaluate the function at the grid points
    eval_points = zeros(ComplexF64, S_pad)

    for n ∈ -div(N_fft,2):div(N_fft,2)-1
        eval_points[n] = G(grid_points[n])
    end

    # Use the FFT to obtain the coefficients of G(u)
    G_of_u.coefficients[:] = project(Sequence(S_pad, [FFTW.fftshift(FFTW.ifft(eval_points.coefficients[begin:end-1]));0]), S).coefficients[:]
end

function G_of_u_vec2vec!(G_of_u::Sequence, u::Sequence, G::Function, N_fft::Int64)
    # INPUTS:
    # 1) u  = {u_n}_{n=-N}^{N} sequence of fourier coefficients
    # 2) g : Rⁿ → Rⁿ analytic function

    # Set ouput to 0
    G_of_u .= 0

    # Extract parameters from inputs
    S = space(component(u,1))
    f = frequency(S)
    dim = div(dimension(space(u)), dimension(S))

    S_pad = Fourier(div(N_fft,2), f)    

    u_pad = project(u, S_pad^dim)

    # Evaluate the fourier series at grid points using the IFFT
    grid_points = zeros(ComplexF64, dim, N_fft)
    for i ∈ 1:dim
        grid_points[i,:] = FFTW.fft(FFTW.ifftshift(component(u_pad,i).coefficients[begin:end-1]))
    end

    # Evaluate the function at the grid points
    eval_points = zeros(ComplexF64, dim, N_fft)

    for j ∈ 1:N_fft
        eval_points[:,j] = G(grid_points[:,j])
    end

    # Use the FFT to obtain the coefficients of G(u)
    fourier_coeffs = zeros(ComplexF64, dim, N_fft)
    for i ∈ 1:dim
        fourier_coeffs[i,:] = FFTW.fftshift(FFTW.ifft(eval_points[i,:]))
    end

    for i ∈ 1:dim
        component(G_of_u,i).coefficients[:] = project(Sequence(S_pad, [fourier_coeffs[i,:];0]), S).coefficients[:]
    end  
end

function G_of_u_vec2pt!(G_of_u::Sequence, u::Sequence, G::Function, N_fft::Int64)
    # INPUTS:
    # 1) u  = {u_n}_{n=-N}^{N} sequence of fourier coefficients
    # 2) g : Rⁿ → R analytic function

    # Set ouput to 0
    G_of_u .= 0

    # Extract parameters from inputs
    S = space(component(u,1))
    f = frequency(S)
    dim = div(dimension(space(u)), dimension(S))

    S_pad = Fourier(div(N_fft,2), f)    

    u_pad = project(u, S_pad^dim)

    # Evaluate the fourier series at grid points using the IFFT
    grid_points = zeros(ComplexF64, dim, N_fft)
    for i ∈ 1:dim
        grid_points[i,:] = FFTW.fft(FFTW.ifftshift(component(u_pad,i).coefficients[begin:end-1]))
    end

    # Evaluate the function at the grid points
    eval_points = zeros(ComplexF64, 1, N_fft)

    for j ∈ 1:N_fft
        eval_points[j] = G(grid_points[:,j])
    end

    # Use the FFT to obtain the coefficients of G(u)
    fourier_coeffs = FFTW.fftshift(FFTW.ifft(eval_points))

    G_of_u.coefficients[:] = project(Sequence(S_pad, [fourier_coeffs;0]), S).coefficients[:]
end

function G_of_u_vec2mat!(G_of_u::LinearOperator, u::Sequence, G::Function, N_fft::Int64)
    # INPUTS:
    # 1) u  = {u_n}_{n=-N}^{N} sequence of fourier coefficients
    # 2) G : Rⁿ → Rⁿ x Rⁿ  analytic function

    # Set ouput to 0
    G_of_u .= 0

    # Extract parameters from inputs
    S = space(component(u,1))
    f = frequency(S)
    dim = div(dimension(space(u)), dimension(S))

    S_pad = Fourier(div(N_fft,2), f)    

    u_pad = project(u, S_pad^dim)

    # Evaluate the fourier series at grid points using the IFFT
    grid_points = zeros(ComplexF64, dim, N_fft)
    for i ∈ 1:dim
        grid_points[i,:] = FFTW.fft(FFTW.ifftshift(component(u_pad,i).coefficients[begin:end-1]))
    end

    # Evaluate the function at the grid points
    eval_points = zeros(ComplexF64, dim, dim, N_fft)

    for j ∈ 1:N_fft
        eval_points[:,:,j] = G(grid_points[:,j])
    end

    # Use the FFT to obtain the coefficients of G(u)
    fourier_coeffs = zeros(ComplexF64, dim, dim, N_fft)
    for i ∈ 1:dim
        for j ∈ 1:dim
            fourier_coeffs[i,j,:] = FFTW.fftshift(FFTW.ifft(eval_points[i,j,:]))
        end
    end

    for i ∈ 1:dim
        for j ∈ 1:dim
            component(G_of_u,i,j).coefficients[:] = project(Multiplication(Sequence(S_pad, [fourier_coeffs[i,j,:];0])), S,S).coefficients[:]
        end
    end     
end

function Newton(u::Sequence, F::Sequence, DF::LinearOperator, tol::Float64 = 1e-12, max_iter::Int64 = 50)
    count = 0;
    F!(F, u, N_fft)
    DF!(DF, u, N_fft)
    while norm(F) > tol && count <= max_iter
        println("Iteration " * string(count) * ", ||F(u)|| = " * string(norm(F)) * ", ||DF\\F|| = " * string(norm(DF \ F)))
        u = u - DF \ F
        F!(F,u,N_fft)
        DF!(DF, u, N_fft)
        count = count + 1
    end
    println("Iteration " * string(count) * ", ||F(u)|| = " * string(norm(F)) * ", ||DF\\F|| = " * string(norm(DF \ F)) * "\nNewton ended after " * string(count) * " iterations needed. \n")
    return u
end

function collection_eval(time_data::Vector{Float64}, u::Sequence)
    num_points = length(time_data)
    space_data = zeros(3, num_points)

    for i ∈ 1:num_points
        space_data[:,i] = real(u(time_data[i]))
    end


    return space_data
end

function kepler_sample(ψ::Vector{Float64},e::Float64, ϕ::Float64, N_fft::Int64)
    # Compting constant c
    f_int(θ,p) = (1+e*cos(θ))^-2
    prob = Integrals.IntegralProblem(f_int, (0.0,2*pi))
    c = (2*pi / Integrals.solve(prob, Integrals.QuadGKJL(); reltol = 1e-14).u)^(1/3)

    # Solving ODE for θ, r, x, y, z
    f_ode(θ,p,t) = c^(-3) * (1+e*cos(θ))^2
    prob = DifferentialEquations.ODEProblem(f_ode, ϕ, (0, 2*pi))
    sol = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5(), reltol = 1e-14, saveat = LinRange(0,2*pi,N_fft+1))
    θ_grid = sol.u[begin:end-1]
    point_grid = zeros(3, N_fft)
    for i ∈ 1:N_fft
        r = c^2 / (1+e*cos(θ_grid[i]))
        point_grid[1,i] = r*cos(θ_grid[i])
        point_grid[2,i] = r*sin(θ_grid[i])
    end
    
    # Building the rotation matrix and applying the rotation
    J₁ = [
        0 0 0
        0 0 -1
        0 1 0
    ]

    J₂ = [
        0 0 -1
        0 0 0
        1 0 0
    ]

    J₃ = [
        0 -1 0
        1 0 0
        0 0 0
    ]

    rot = exp(Ψ[1]*J₁ + Ψ[2]*J₂ + Ψ[3]*J₃)

    return rot * point_grid, sol.t[begin:end-1]
end

function generators()
    A = [
        1 0 0
        0 -1 0
        0 0 -1
    ]

    B = [
        0 1 0
        0 0 1
        1 0 0
    ]

    generators = zeros(3,3, 11)
    generators[:,:,1] = A
    generators[:,:,2] = B
    generators[:,:,3] = B^2
    generators[:,:,4] = A*B
    generators[:,:,5] = A*B^2
    generators[:,:,6] = B*A
    generators[:,:,7] = B^2*A
    generators[:,:,8] = A*B*A
    generators[:,:,9] = A*B^2*A
    generators[:,:,10] = B*A*B
    generators[:,:,11] = B^2*A*B

    return generators, 11
end