function F!(F::Sequence, u::Sequence, ε::Float64, N_fft::Int64)
    # Extract important space data
    S = space(component(u,1))
    G = zeros(ComplexF64, S^3)   
    N = order(S)

    # Initializing the output to be 0
    F .= 0

    # Linear part
    D = zeros(S^3,S^3)
    for i = 1:3
        for k in -N:N
            component(D,i,i)[k,k] = -k^2
        end
    end

    # Nonlinear part
    G_of_u_vec2vec!(G, u, g, N_fft)

    if ε == 0
        F[:] = (D*u + G)[:]
    else
        # Pertubation part
        H = zeros(ComplexF64, S^3) 
        G_of_u_vec2vec!(H, u, h, N_fft)

        F[:] = (D*u + G + ε*H)[:] 
    end
end

function DF!(DF::LinearOperator, u::Sequence, ε::Float64, N_fft::Int64)
    # Extract important space data
    S = space(component(u,1))
    DG = zeros(ComplexF64, S^3, S^3)   
    N = order(S)

    # Initializing the output to be 0
    DF .= 0

    # Linear part
    D = zeros(S^3,S^3)
    for i = 1:3
        for k in -N:N
            component(D,i,i)[k,k] = -k^2
        end
    end

    # Nonlinear part
    G_of_u_vec2mat!(DG, u, Dg, N_fft)

    if ε == 0
        DF.coefficients[:] = (D + DG).coefficients[:]
    else
        # Pertubation part
        DH = zeros(ComplexF64, S^3, S^3) 
        G_of_u_vec2mat!(DH, u, Dh, N_fft)

        DF.coefficients[:] = (D + DG + ε*DH).coefficients[:]
    end
end

function g(u)
    return u/(u[1]^2 + u[2]^2 + u[3]^2)^(3/2)
end

function Dg(u)
    n = (u[1]^2 + u[2]^2 + u[3]^2)^(1/2)
    return 1/n^3 * I -3*u*transpose(u)/n^5
end

function h(u)
    gens, num_gen = generators()
    sum = zeros(3,1)    

    for i ∈ 1:num_gen
        A = I - gens[:,:,i]
        Au = A*u
        sum = sum + Au/(Au[1]^2 + Au[2]^2 + Au[3]^2)^(3/2)
    end

    return sum
end

function Dh(u)
    gens, num_gen = generators()
    sum = zeros(3,3)    

    for i ∈ 1:num_gen
        A = I - gens[:,:,i]
        Au = A*u
        n = (Au[1]^2 + Au[2]^2 + Au[3]^2)^(1/2)
        sum = sum + A/n^3  - 3*Au*transpose(Au)*A/n^5
    end

    return sum
end

function A_1(ψ, e, N_fft)
    sample_points, sample_time = kepler_sample(ψ, e, 0.0, N_fft)

    function f_A(u)
        gens, num_gen = generators()
        sum = 0    

        for i ∈ 1:num_gen
            A = I - gens[:,:,i]
            Au = A*u

            # if norm(Au) < 1e-6
            #     println("Particles are too close\nProblem occured using u =" * string(u) * "and generator g =" * string(gens[:,:,i]))
            # end
            sum = sum + 1/((Au[1]^2 + Au[2]^2 + Au[3]^2)^(1/2))
            # println(sum)
        end

        return sum
    end
    f_points = zeros(1, N_fft)
    for i ∈ 1:N_fft
        f_points[i] = f_A(sample_points[:,i])
    end

    # display(f_points)

    return 2*pi/N_fft * FFTW.fft(f_points)[1]
end

function A_1_theta(ψ, e, N_fft)
    # Compting constant c
    f_int(θ,p) = (1+e*cos(θ))^-2
    prob = Integrals.IntegralProblem(f_int, (0.0,2*pi))
    c = (2*pi / Integrals.solve(prob, Integrals.QuadGKJL(); reltol = 1e-14).u)^(1/3)

    function f_A(θ)
        gens, num_gen = generators()
        sum = 0

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

        rot = exp(ψ[3]*J₃) * exp(ψ[2]*J₂) * exp(ψ[1]*J₁)

        u = rot * [cos(θ); sin(θ); 0]    

        for i ∈ 1:num_gen
            A = I - gens[:,:,i]
            Au = A*u

            sum = sum + 1/((Au[1]^2 + Au[2]^2 + Au[3]^2)^(1/2))
        end

        return sum
    end
    f_points = zeros(1, N_fft)
    for i ∈ 1:N_fft
        θ = 2*pi*(i-1) / N_fft
        f_points[i] = 1/(1+e*cos(θ)) * f_A(θ)
    end

    return c*2*pi/N_fft * FFTW.fft(f_points)[1]
end