include("integrate.jl")

function generators_12()
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

    gen = zeros(Int64, 3, 3, 12)
    gen[:, :, 1] = A^0
    gen[:, :, 2] = B
    gen[:, :, 3] = B^2
    gen[:, :, 4] = A
    gen[:, :, 5] = A * B
    gen[:, :, 6] = A * B^2
    gen[:, :, 7] = B * A
    gen[:, :, 8] = B^2 * A
    gen[:, :, 9] = A * B * A
    gen[:, :, 10] = B * A * B
    gen[:, :, 11] = B^2 * A * B
    gen[:, :, 12] = A * B^2 * A * B

    return gen, 12
end

function generators_24()
    A = [
        0 -1 0
        1 0 0
        0 0 1
    ]

    B = [
        -1 -0 0
        0 0 1
        1 0 0
    ]

    gen = zeros(Int64, 3, 3, 24)

    # for i ∈ 0:3
    #     gen[:, :, 6*i + 1] = A^i * I
    #     gen[:, :, 6*i + 2] = A^i * (B)
    #     gen[:, :, 6*i + 3] = A^i * (B * A)
    #     gen[:, :, 6*i + 4] = A^i * (B * A * B)
    #     gen[:, :, 6*i + 5] = A^i * (B * A * B * A)
    #     gen[:, :, 6*i + 6] = A^i * (B * A * B * A * B)
    # end

    # gen[:, :, 1]  = A^0
    # gen[:, :, 2]  = A
    # gen[:, :, 3]  = B
    # gen[:, :, 4]  = A^2
    # gen[:, :, 5]  = B*A
    # gen[:, :, 6]  = A*B
    # gen[:, :, 7]  = A^3
    # gen[:, :, 8]  = B*A^2
    # gen[:, :, 9]  = A*B*A
    # gen[:, :, 10] = A^2*B
    # gen[:, :, 11] = B*A*B
    # gen[:, :, 12] = B*A^3
    # gen[:, :, 13] = A*B*A^2
    # gen[:, :, 14] = A^2*B*A
    # gen[:, :, 15] = B*A*B*A
    # gen[:, :, 16] = B*A^2*B
    # gen[:, :, 17] = A*B*A^3
    # gen[:, :, 18] = A^2*B*A^2
    # gen[:, :, 19] = B*A*B*A^2
    # gen[:, :, 20] = B*A^2*B*A
    # gen[:, :, 21] = A*B*A^2*B
    # gen[:, :, 22] = A^2*B*A^3
    # gen[:, :, 23] = B*A*B*A^3
    # gen[:, :, 24] = B*A^2*B*A^2

    gen[:,:,1]  = [ 1  0  0;  0  1  0;  0  0  1]

    gen[:,:,2]  = [ 1  0  0;  0 -1  0;  0  0 -1]
    gen[:,:,3]  = [-1  0  0;  0  1  0;  0  0 -1]
    gen[:,:,4]  = [-1  0  0;  0 -1  0;  0  0  1]

    gen[:,:,5]  = [0  1  0;  1  0  0;  0  0 -1]
    gen[:,:,6]  = [0  1  0; -1  0  0;  0  0  1]
    gen[:,:,7]  = [0 -1  0;  1  0  0;  0  0  1]
    gen[:,:,8]  = [0 -1  0; -1  0  0;  0  0 -1]

    gen[:,:,9]  = [0  0  1;  0  1  0; -1  0  0]
    gen[:,:,10] = [0  0  1;  0 -1  0;  1  0  0]
    gen[:,:,11] = [0  0 -1;  0  1  0;  1  0  0]
    gen[:,:,12] = [0  0 -1;  0 -1  0; -1  0  0]

    gen[:,:,13] = [0  1  0;  0  0  1;  1  0  0]
    gen[:,:,14] = [0  1  0;  0  0 -1; -1  0  0]
    gen[:,:,15] = [0 -1  0;  0  0  1; -1  0  0]
    gen[:,:,16] = [0 -1  0;  0  0 -1;  1  0  0]

    gen[:,:,17] = [0  0  1;  1  0  0;  0  1  0]
    gen[:,:,18] = [0  0  1; -1  0  0;  0 -1  0]
    gen[:,:,19] = [0  0 -1;  1  0  0;  0 -1  0]
    gen[:,:,20] = [0  0 -1; -1  0  0;  0  1  0]

    gen[:,:,21] = [ 1  0  0; 0  0 -1; 0  1  0]
    gen[:,:,22] = [ 1  0  0; 0  0  1; 0 -1  0]
    gen[:,:,23] = [-1  0  0; 0  0  1; 0  1  0]
    gen[:,:,24] = [-1  0  0; 0  0 -1; 0 -1  0]

    return gen, 24
end

function generators_60()
    Φ = 1/2 * (1 + sqrt(5))

    A = 1/2 * [
        Φ-1  Φ   1
        Φ    -1  Φ-1
        1    Φ-1 -Φ
    ]

    B = [
         0 0 1
         1 0 0
         0 1 0
    ]

    gen = zeros(3, 3, 60)

    # gen[:, :, 1] = A^0
    # gen[:, :, 2] = A
    # gen[:, :, 3] = B
    # gen[:, :, 4] = B * A
    # gen[:, :, 5] = A * B
    # gen[:, :, 6] = B^2

    # gen[:, :, 7] = A * B * A
    # gen[:, :, 8] = B^2 * A
    # gen[:, :, 9] = B * A * B
    # gen[:, :, 10] = A * B^2
    # gen[:, :, 11] = B * A * B * A
    # gen[:, :, 12] = A * B^2 * A
    # gen[:, :, 13] = A * B * A * B
    # gen[:, :, 14] = B^2 * A * B
    # gen[:, :, 15] = B * A * B^2

    # gen[:, :, 16] = A * B * A * B * A
    # gen[:, :, 17] = B^2 * A * B * A
    # gen[:, :, 18] = B * A * B^2 * A
    # gen[:, :, 19] = B * A * B * A * B
    # gen[:, :, 20] = A * B^2 * A * B
    # gen[:, :, 21] = A * B * A * B^2
    # gen[:, :, 22] = B^2 * A * B^2

    # gen[:, :, 23] = B * A * B * A * B * A
    # gen[:, :, 24] = A * B^2 * A * B * A
    # gen[:, :, 25] = A * B * A * B^2 * A
    # gen[:, :, 26] = B^2 * A * B^2 * A
    # gen[:, :, 27] = B^2 * A * B * A * B
    # gen[:, :, 28] = B * A * B^2 * A * B
    # gen[:, :, 29] = B * A * B * A * B^2

    # gen[:, :, 30] = B^2 * A * B * A * B * A
    # gen[:, :, 31] = B * A * B^2 * A * B * A
    # gen[:, :, 32] = B * A * B * A * B^2 * A
    # gen[:, :, 33] = A * B^2 * A * B * A * B
    # gen[:, :, 34] = A * B * A * B^2 * A * B
    # gen[:, :, 35] = B^2 * A * B^2 * A * B
    # gen[:, :, 36] = B^2 * A * B * A * B^2

    # gen[:, :, 37] = A * B^2 * A * B * A * B * A
    # gen[:, :, 38] = A * B * A * B^2 * A * B * A
    # gen[:, :, 39] = B^2 * A * B^2 * A * B * A
    # gen[:, :, 40] = B^2 * A * B * A * B^2 * A
    # gen[:, :, 41] = B * A * B^2 * A * B * A * B
    # gen[:, :, 42] = B * A * B * A * B^2 * A * B
    # gen[:, :, 43] = A * B^2 * A * B * A * B^2

    # gen[:, :, 44] = B * A * B^2 * A * B * A * B * A
    # gen[:, :, 45] = B * A * B * A * B^2 * A * B * A
    # gen[:, :, 46] = A * B^2 * A * B * A * B^2 * A
    # gen[:, :, 47] = A * B * A * B^2 * A * B * A * B
    # gen[:, :, 48] = B^2 * A * B^2 * A * B * A * B
    # gen[:, :, 49] = B^2 * A * B * A * B^2 * A * B
    # gen[:, :, 50] = B * A * B^2 * A * B * A * B^2

    # gen[:, :, 51] = A * B * A * B^2 * A * B * A * B * A
    # gen[:, :, 52] = B^2 * A * B^2 * A * B * A * B * A
    # gen[:, :, 53] = B^2 * A * B * A * B^2 * A * B * A
    # gen[:, :, 54] = B * A * B^2 * A * B * A * B^2 * A
    # gen[:, :, 55] = A * B^2 * A * B * A * B^2 * A * B
    # gen[:, :, 56] = A * B * A * B^2 * A * B * A * B^2

    # gen[:, :, 57] = A * B^2 * A * B * A * B^2 * A * B * A
    # gen[:, :, 58] = A * B * A * B^2 * A * B * A * B^2 * A
    # gen[:, :, 59] = B * A * B^2 * A * B * A * B^2 * A * B
    # gen[:, :, 60] = B * A * B^2 * A * B * A * B^2 * A * B * A

    gen[:,:,1] = A^0
    gen[:,:,2] = A
    gen[:,:,3] = B
    gen[:,:,4] = B*A
    gen[:,:,5] = A*B
    gen[:,:,6] = B^2
    gen[:,:,7] = A*B*A
    gen[:,:,8] = B^2*A
    gen[:,:,9] = B*A*B
    gen[:,:,10] = A*B^2
    gen[:,:,11] = B*A*B*A
    gen[:,:,12] = A*B^2*A
    gen[:,:,13] = A*B*A*B
    gen[:,:,14] = B^2*A*B
    gen[:,:,15] = B*A*B^2
    gen[:,:,16] = A*B*A*B*A
    gen[:,:,17] = B^2*A*B*A
    gen[:,:,18] = B*A*B^2*A
    gen[:,:,19] = B*A*B*A*B
    gen[:,:,20] = A*B^2*A*B
    gen[:,:,21] = A*B*A*B^2
    gen[:,:,22] = B^2*A*B^2
    gen[:,:,23] = B*A*B*A*B*A
    gen[:,:,24] = A*B^2*A*B*A
    gen[:,:,25] = A*B*A*B^2*A
    gen[:,:,26] = B^2*A*B^2*A
    gen[:,:,27] = B^2*A*B*A*B
    gen[:,:,28] = B*A*B^2*A*B
    gen[:,:,29] = B*A*B*A*B^2
    gen[:,:,30] = B^2*A*B*A*B*A
    gen[:,:,31] = B*A*B^2*A*B*A
    gen[:,:,32] = B*A*B*A*B^2*A
    gen[:,:,33] = A*B^2*A*B*A*B
    gen[:,:,34] = A*B*A*B^2*A*B
    gen[:,:,35] = B^2*A*B^2*A*B
    gen[:,:,36] = B^2*A*B*A*B^2
    gen[:,:,37] = A*B^2*A*B*A*B*A
    gen[:,:,38] = A*B*A*B^2*A*B*A
    gen[:,:,39] = B^2*A*B^2*A*B*A
    gen[:,:,40] = B^2*A*B*A*B^2*A
    gen[:,:,41] = B*A*B^2*A*B*A*B
    gen[:,:,42] = B*A*B*A*B^2*A*B
    gen[:,:,43] = A*B^2*A*B*A*B^2
    gen[:,:,44] = B*A*B^2*A*B*A*B*A
    gen[:,:,45] = B*A*B*A*B^2*A*B*A
    gen[:,:,46] = A*B^2*A*B*A*B^2*A
    gen[:,:,47] = A*B*A*B^2*A*B*A*B
    gen[:,:,48] = B^2*A*B^2*A*B*A*B
    gen[:,:,49] = B^2*A*B*A*B^2*A*B
    gen[:,:,50] = B*A*B^2*A*B*A*B^2
    gen[:,:,51] = A * B * A * B^2 * A * B * A * B * A         # ABAB²ABABA
    gen[:,:,52] = B^2 * A * B^2 * A * B * A * B * A           # B²AB²ABABA
    gen[:,:,53] = B^2 * A * B * A * B^2 * A * B * A           # B²ABAB²ABA
    gen[:,:,54] = B * A * B^2 * A * B * A * B^2 * A           # BAB²ABAB²A
    gen[:,:,55] = A * B^2 * A * B * A * B^2 * A * B           # AB²ABAB²AB
    gen[:,:,56] = A * B * A * B^2 * A * B * A * B^2           # ABAB²ABAB²
    gen[:,:,57] = A * B^2 * A * B * A * B^2 * A * B * A       # AB²ABAB²ABA
    gen[:,:,58] = A * B * A * B^2 * A * B * A * B^2 * A       # ABAB²ABAB²A
    gen[:,:,59] = B * A * B^2 * A * B * A * B^2 * A * B       # BAB²ABAB²AB
    gen[:,:,60] = B * A * B^2 * A * B * A * B^2 * A * B * A   # BAB²ABAB²ABA

    return gen, 60
end

function rotation_matrices()
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

    return (J₁, J₂, J₃)
end

function kepler_sample(ψ::Vector{Float64}, e::Float64, ϕ::Float64, num_points::Int64)
    # Computing constant c
    f_c = θ -> (1 + e * cos(θ))^(-2)
    c = (2 * ((π) / fft_integrate_float(f_c, 2 //1, 2^14)))^(1 / 3)

    # Solving ODE for θ, which lets us compute r, x, y, z
    f_ode(θ, p, t) = c^(-3) * (1 + e * cos(θ))^2
    prob = ODEProblem(f_ode, ϕ, (0, 2 * pi))
    sol = solve(prob, DifferentialEquations.Tsit5(), reltol=1e-14, abstol=1e-14, saveat=LinRange(0, 2 * pi, num_points + 1))
    θ_grid = sol.u[begin:end-1]
    point_grid = zeros(3, num_points)
    for i ∈ 1:num_points
        r = c^2 / (1 + e * cos(θ_grid[i]))
        point_grid[1, i] = r * cos(θ_grid[i])
        point_grid[2, i] = r * sin(θ_grid[i])
    end

    # Defining rotation matrices
    J₁, J₂, J₃ = rotation_matrices()

    𝒥₁ = [
        1 0 0
        0 cos(ψ[1]) -sin(ψ[1])
        0 sin(ψ[1]) cos(ψ[1])
    ]
    𝒥₂ = [
        cos(ψ[2]) 0 -sin(ψ[2])
        0 1 0
        sin(ψ[2]) 0 cos(ψ[2])
    ]
    𝒥₃ = [
        cos(ψ[3]) -sin(ψ[3]) 0
        sin(ψ[3]) cos(ψ[3]) 0
        0 0 1
    ]

    rot = 𝒥₃ * 𝒥₂ * 𝒥₁

    return rot * point_grid, sol.t[begin:end-1]
end

function collection_eval(time_data::Vector{Float64}, u::Sequence)
    num_points = length(time_data)
    space_data = zeros(3, num_points)

    for i ∈ 1:num_points
        space_data[:, i] = real(u(time_data[i]))
    end

    return space_data
end

