include("fourier_computations.jl")

function F!(F::Sequence, u::Sequence, ε::Float64, N_fft::Int64)
    # Extract important space data
    S = space(block(u, 1))
    G = zeros(ComplexF64, S^3)
    N = order(S)

    # Initializing the output to be 0
    F.coefficients .= 0

    # Linear part
    D = zeros(S^3, S^3)
    for i = 1:3
        for k in -N:N
            block(D, i, i)[k, k] = -k^2
        end
    end

    # Nonlinear part
    G_of_u_vec2vec!(G, u, g, N_fft)

    if ε == 0
        F[:] = (D*u+G)[:]
    else
        # Pertubation part
        H = zeros(ComplexF64, S^3)
        G_of_u_vec2vec!(H, u, h, N_fft)

        F[:] = (D*u+G+ε*H)[:]
    end
end

function DF!(DF::LinearOperator, u::Sequence, ε::Float64, N_fft::Int64)
    # Extract important space data
    S = space(block(u, 1))
    DG = zeros(ComplexF64, S^3, S^3)
    N = order(S)

    # Initializing the output to be 0
    DF.coefficients .= 0

    # Linear part
    D = zeros(S^3, S^3)
    for i = 1:3
        for k in -N:N
            block(D, i, i)[k, k] = -k^2
        end
    end

    # Nonlinear part
    G_of_u_vec2mat!(DG, u, Dg, N_fft)

    if ε == 0
        DF.coefficients[:] = (D+DG).coefficients[:]
    else
        # Pertubation part
        DH = zeros(ComplexF64, S^3, S^3)
        G_of_u_vec2mat!(DH, u, Dh, N_fft)

        DF.coefficients[:] = (D+DG+ε*DH).coefficients[:]
    end
end

function g(u)
    return u / (u[1]^2 + u[2]^2 + u[3]^2)^(3 / 2)
end

function Dg(u)
    n = (u[1]^2 + u[2]^2 + u[3]^2)^(1 / 2)
    return 1 / n^3 * I - 3 * u * transpose(u) / n^5
end

function h(u)
    gens, num_gen = generators()
    sum = zeros(3, 1)

    for i ∈ 1:num_gen
        A = I - gens[:, :, i]
        Au = A * u
        sum = sum + Au / (Au[1]^2 + Au[2]^2 + Au[3]^2)^(3 / 2)
    end

    return sum
end

function Dh(u)
    gens, num_gen = generators()
    sum = zeros(3, 3)

    for i ∈ 1:num_gen
        A = I - gens[:, :, i]
        Au = A * u
        n = (Au[1]^2 + Au[2]^2 + Au[3]^2)^(1 / 2)
        sum = sum + A / n^3 - 3 * Au * transpose(Au) * A / n^5
    end

    return sum
end