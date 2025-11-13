function F!(F::Sequence, u::Sequence, N_fft::Int64)
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

    F[:] = (D*u + G)[:] 
end

function DF!(DF::LinearOperator, u::Sequence, N_fft::Int64)
    # Extract important space data
    S = space(component(u,1))
    DG = zeros(ComplexF64, S^3,S^3)    
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

    DF.coefficients[:] = (D + DG).coefficients[:]  
end

function g(u)
    return u/((u[1]^2 + u[2]^2 + u[3]^2)^(3/2))
end

function Dg(u)
    n = (u[1]^2 + u[2]^2 + u[3]^2)^(1/2)
    return 1/(n^3) * I -3*u*transpose(u)/(n^5)
end

function h(u)
    gens, num_gen = generators()
    sum = zeros(3,1)    

    for i ∈ 1:num_gen
        A = I - gens[:,:,i]
        Au = A*u
        sum = sum + (Au)/((Au[1]^2 + Au[2]^2 + Au[3]^2)^(3/2))
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

            if norm(Au) < 1e-6
                println("Particles are too close\nProblem occured using u =" * string(u) * "and generator g =" * string(gens[:,:,i]))
            end
            sum = sum + 1/((Au[1]^2 + Au[2]^2 + Au[3]^2)^(1/2))
            println(sum)
        end

        return sum
    end
    f_points = zeros(1, N_fft)
    for i ∈ 1:N_fft
        f_points[i] = f_A(sample_points[:,i])
    end

    display(f_points)

    return 2*pi/N_fft * FFTW.fft(f_points)[1]
end