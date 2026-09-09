include("integrate.jl")
include("helpers.jl")

function f_A(u::Vector{Float64})
    gens, num_gen = generators()
    tmp = 0

    println("u = " * string(u))

    for i ∈ 1:num_gen
        A = I - gens[:, :, i]
        Au = A * u

        println("u - Lu = " * string(Au))

        tmp = tmp + 1 / ((Au[1]^2 + Au[2]^2 + Au[3]^2)^(1 / 2))
    end
    println()

    return tmp
end

function f_A(u::Vector{Complex{Interval{Float64}}})
    gens, num_gen = generators()
    tmp = interval(0)

    for i ∈ 1:num_gen
        A = interval.(I - gens[:, :, i])
        Au = A * u

        tmp = tmp + interval(1) / ((Au[1]^2 + Au[2]^2 + Au[3]^2)^(interval(1) / interval(2)))
    end

    return tmp
end

function f_DA(u::Vector{Float64}, du::Vector{Float64})
    gens, num_gen = generators()
    tmp = 0

    for i ∈ 1:num_gen
        A = I - gens[:, :, i]
        Au = A * u
        Adu = A * du

        tmp = tmp - 1 / ((Au[1] .^ 2 + Au[2]^2 + Au[3]^2)^(3 / 2)) * sum(Au .* Adu)
    end

    return tmp
end

function f_DA(u::Vector{Complex{Interval{Float64}}}, du::Vector{Complex{Interval{Float64}}})
    gens, num_gen = generators()
    tmp = interval(0)

    for i ∈ 1:num_gen
        A = interval.(I - gens[:, :, i])
        Au = A * u
        Adu = A * du

        tmp = tmp - exact(1) / ((Au[1] .^ 2 + Au[2]^2 + Au[3]^2)^(interval(3) / interval(2))) * sum(Au .* Adu)
    end

    return tmp
end

function f_HA(u::Vector{Float64}, du_1::Vector{Float64}, du_2::Vector{Float64}, du_12::Vector{Float64})
    gens, num_gen = generators()
    tmp = 0

    for i ∈ 1:num_gen
        A = I - gens[:, :, i]
        Au = A * u
        Adu_1 = A * du_1
        Adu_2 = A * du_2
        Adu_12 = A * du_12

        tmp = tmp - 1 / ((Au[1]^2 + Au[2]^2 + Au[3]^2)^(3 / 2)) * sum(Adu_2 .* Adu_1)
        tmp = tmp - 1 / ((Au[1]^2 + Au[2]^2 + Au[3]^2)^(3 / 2)) * sum(Au .* Adu_12)
        tmp = tmp + 3 / ((Au[1]^2 + Au[2]^2 + Au[3]^2)^(5 / 2)) * sum(Au .* Adu_1) * sum(Au .* Adu_2)
    end

    return tmp
end

function f_HA(u::Vector{Complex{Interval{Float64}}}, du_1::Vector{Complex{Interval{Float64}}}, du_2::Vector{Complex{Interval{Float64}}}, du_12::Vector{Complex{Interval{Float64}}})
    gens, num_gen = generators()
    tmp = interval(0)

    for i ∈ 1:num_gen
        A = interval.(I - gens[:, :, i])
        Au = A * u
        Adu_1 = A * du_1
        Adu_2 = A * du_2
        Adu_12 = A * du_12

        tmp = tmp - interval(1) / ((Au[1]^2 + Au[2]^2 + Au[3]^2)^interval(3) / interval(2)) * sum(Adu_2 .* Adu_1)
        tmp = tmp - interval(1) / ((Au[1]^2 + Au[2]^2 + Au[3]^2)^interval(3) / interval(2)) * sum(Au .* Adu_12)
        tmp = tmp + interval(3) / ((Au[1]^2 + Au[2]^2 + Au[3]^2)^interval(5) / interval(2)) * sum(Au .* Adu_1) * sum(Au .* Adu_2)
    end

    return tmp
end

function A_1(X::Vector{Float64}, N_fft::Int64)
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    # Computing constant c
    f_c = θ -> (1 + e * cos(θ))^-2
    c = (2 * (pi / fft_integrate(f_c, 2 * pi, N_fft)))^(1 / 3)

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

    u_base = θ -> [cos(θ); sin(θ); 0]
    u = θ -> exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁) * u_base(θ)

    return c * fft_integrate(θ -> f_A(u(θ)) * 1 / (1 + e * cos(θ)), 2 * π, N_fft)
end

function A_1(X::Vector{Interval{Float64}}, N_fft::Int64)
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    # Computing constant c
    @exact f_c = θ -> (1 + e * cos(θ))^-2
    c = (interval(2) * (pi / fft_integrate(f_c, interval(2) * pi, N_fft)))^(interval(1) / interval(3))

    J₁ = interval.([
        0 0 0
        0 0 -1
        0 1 0
    ])

    J₂ = interval.([
        0 0 -1
        0 0 0
        1 0 0
    ])

    J₃ = interval.([
        0 -1 0
        1 0 0
        0 0 0
    ])

    u_base = θ -> [cos(θ); sin(θ); interval(0)]
    u = θ -> exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁) * u_base(θ)

    return c * fft_integrate(θ -> f_A(u(θ)) * interval(1) / (interval(1) + e * cos(θ)), interval(2) * π, N_fft; rigorous=true)
end


function DA_1(X::Vector{Float64}, N_fft::Int64)
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    derivative = zeros(4, 1)

    # Computing constant c
    f_c = θ -> (1 + e * cos(θ))^-2
    c = (2 * pi / fft_integrate(f_c, 2 * pi, N_fft))^(1 / 3)

    # Defining rotation matrices
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

    u_base = θ -> [cos(θ); sin(θ); 0]

    u = θ -> (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)

    du_1 = θ -> (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * J₁ * exp(ψ[1] * J₁)) * u_base(θ)
    du_2 = θ -> (exp(ψ[3] * J₃) * J₂ * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)
    du_3 = θ -> (J₃ * exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)

    # Computing derivatives with respect to ψ  

    # ψ₁
    derivative[1] = c * fft_integrate(θ -> f_DA(u(θ), du_1(θ)) * 1 / (1 + e * cos(θ)), 2 * π, N_fft)
    # ψ₂
    derivative[2] = c * fft_integrate(θ -> f_DA(u(θ), du_2(θ)) * 1 / (1 + e * cos(θ)), 2 * π, N_fft)
    # ψ₃
    derivative[3] = c * fft_integrate(θ -> f_DA(u(θ), du_3(θ)) * 1 / (1 + e * cos(θ)), 2 * π, N_fft)

    # Computing derivative with respect to e 

    dc = c^4 / (3 * pi) * fft_integrate(θ -> cos(θ) / ((1 + e * cos(θ))^3), 2 * π, N_fft)
    derivative[4] = fft_integrate(θ -> (dc / (1 + e * cos(θ)) - c * cos(θ) / ((1 + e * cos(θ))^2)) * f_A(u(θ)), 2 * π, N_fft)

    return derivative
end

function DA_1(X::Vector{Interval{Float64}}, N_fft::Int64)
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    derivative = zeros(Interval{Float64}, 4, 1)

    # Computing constant c
    @exact f_c = θ -> (1 + e * cos(θ))^-2
    c = (interval(2) * pi / fft_integrate(f_c, interval(2) * pi, N_fft))^(interval(1) / interval(3))

    # Defining rotation matrices
    J₁ = interval.([
        0 0 0
        0 0 -1
        0 1 0
    ])

    J₂ = interval.([
        0 0 -1
        0 0 0
        1 0 0
    ])

    J₃ = interval.([
        0 -1 0
        1 0 0
        0 0 0
    ])

    u_base = θ -> [cos(θ); sin(θ); interval(0)]

    u = θ -> (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)

    du_1 = θ -> (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * J₁ * exp(ψ[1] * J₁)) * u_base(θ)
    du_2 = θ -> (exp(ψ[3] * J₃) * J₂ * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)
    du_3 = θ -> (J₃ * exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)

    # Computing derivatives with respect to ψ  

    # ψ₁
    derivative[1] = c * fft_integrate(θ -> f_DA(u(θ), du_1(θ)) * interval(1) / (interval(1) + e * cos(θ)), interval(2) * π, N_fft, rigorous=true)
    # ψ₂
    derivative[2] = c * fft_integrate(θ -> f_DA(u(θ), du_2(θ)) * interval(1) / (interval(1) + e * cos(θ)), interval(2) * π, N_fft, rigorous=true)
    # ψ₃
    derivative[3] = c * fft_integrate(θ -> f_DA(u(θ), du_3(θ)) * interval(1) / (interval(1) + e * cos(θ)), interval(2) * π, N_fft, rigorous=true)

    # Computing derivative with respect to e 

    dc = c^4 / (interval(3) * π) * fft_integrate(θ -> cos(θ) / ((interval(1) + e * cos(θ))^3), interval(2) * π, N_fft, rigorous=true)
    derivative[4] = fft_integrate(θ -> (dc / (interval(1) + e * cos(θ)) - c * cos(θ) / ((interval(1) + e * cos(θ))^2)) * f_A(u(θ)), interval(2) * π, N_fft)

    return derivative
end

function HA_1(X::Vector{Float64}, N_fft::Int64)
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    derivative = zeros(4, 4)

    # Computing constant c
    f_c = θ -> (1 + e * cos(θ))^-2
    c = (2 * pi / fft_integrate(f_c, 2 * pi, N_fft))^(1 / 3)

    # Defining rotation matrices
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

    u_base = θ -> [cos(θ); sin(θ); 0]

    u = θ -> (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)

    du_1 = θ -> (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * J₁ * exp(ψ[1] * J₁)) * u_base(θ)
    du_2 = θ -> (exp(ψ[3] * J₃) * J₂ * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)
    du_3 = θ -> (J₃ * exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)

    du_11 = θ -> (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * J₁^2 * exp(ψ[1] * J₁)) * u_base(θ)
    du_12 = θ -> (exp(ψ[3] * J₃) * J₂ * exp(ψ[2] * J₂) * J₁ * exp(ψ[1] * J₁)) * u_base(θ)
    du_13 = θ -> (J₃ * exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * J₁ * exp(ψ[1] * J₁)) * u_base(θ)

    du_22 = θ -> (exp(ψ[3] * J₃) * J₂^2 * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)
    du_23 = θ -> (J₃ * exp(ψ[3] * J₃) * J₂ * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)

    du_33 = θ -> (J₃^2 * exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)

    # Computing derivatives with respect to ψ  

    derivative[1, 1] = c * fft_integrate(θ -> f_HA(u(θ), du_1(θ), du_1(θ), du_11(θ)) * 1 / (1 + e * cos(θ)), 2 * π, N_fft)
    derivative[1, 2] = c * fft_integrate(θ -> f_HA(u(θ), du_1(θ), du_2(θ), du_12(θ)) * 1 / (1 + e * cos(θ)), 2 * π, N_fft)
    derivative[1, 3] = c * fft_integrate(θ -> f_HA(u(θ), du_1(θ), du_3(θ), du_13(θ)) * 1 / (1 + e * cos(θ)), 2 * π, N_fft)

    derivative[2, 1] = derivative[1, 2]
    derivative[2, 2] = c * fft_integrate(θ -> f_HA(u(θ), du_2(θ), du_2(θ), du_22(θ)) * 1 / (1 + e * cos(θ)), 2 * π, N_fft)
    derivative[2, 3] = c * fft_integrate(θ -> f_HA(u(θ), du_2(θ), du_3(θ), du_23(θ)) * 1 / (1 + e * cos(θ)), 2 * π, N_fft)

    derivative[3, 1] = derivative[1, 3]
    derivative[3, 2] = derivative[2, 3]
    derivative[3, 3] = c * fft_integrate(θ -> f_HA(u(θ), du_3(θ), du_3(θ), du_33(θ)) * 1 / (1 + e * cos(θ)), 2 * π, N_fft)

    # Computing derivatives with respect to e 

    dc = c^4 / (3 * pi) * fft_integrate(θ -> cos(θ) / ((1 + e * cos(θ))^3), 2 * π, N_fft)
    d2c = 4 * c^3 * dc / (3 * pi) * fft_integrate(θ -> cos(θ) / ((1 + e * cos(θ))^3), 2 * π, N_fft) -
          c^4 / pi * fft_integrate(θ -> cos(θ)^2 / ((1 + e * cos(θ))^4), 2 * π, N_fft)

    de_f = θ -> dc / (1 + e * cos(θ)) - c * cos(θ) / ((1 + e * cos(θ))^2)
    de_2_f = θ -> d2c / (1 + e * cos(θ)) -
                  2 * dc * cos(θ) / ((1 + e * cos(θ))^2) +
                  2 * c * (cos(θ)^2) / ((1 + e * cos(θ))^3)


    derivative[4, 4] = fft_integrate(θ -> f_A(u(θ)) * de_2_f(θ), 2 * π, N_fft)

    # Computing mixed derivatives
    derivative[1, 4] = fft_integrate(θ -> f_DA(u(θ), du_1(θ)) * de_f(θ), 2 * π, N_fft)
    derivative[4, 1] = derivative[1, 4]

    derivative[2, 4] = fft_integrate(θ -> f_DA(u(θ), du_2(θ)) * de_f(θ), 2 * π, N_fft)
    derivative[4, 2] = derivative[2, 4]

    derivative[3, 4] = fft_integrate(θ -> f_DA(u(θ), du_3(θ)) * de_f(θ), 2 * π, N_fft)
    derivative[4, 3] = derivative[3, 4]

    return derivative
end

function HA_1(X::Vector{Interval{Float64}}, N_fft::Int64)
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    derivative = zeros(Interval{Float64}, 4, 4)

    # Computing constant c
    @exact f_c = θ -> (1 + e * cos(θ))^-2
    c = (interval(2) * pi / fft_integrate(f_c, interval(2) * pi, N_fft))^(interval(1) / interval(3))

    # Defining rotation matrices
    J₁ = interval.([
        0 0 0
        0 0 -1
        0 1 0
    ])

    J₂ = interval.([
        0 0 -1
        0 0 0
        1 0 0
    ])

    J₃ = interval.([
        0 -1 0
        1 0 0
        0 0 0
    ])

    u_base = θ -> [cos(θ); sin(θ); interval(0)]

    u = θ -> (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)

    du_1 = θ -> (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * J₁ * exp(ψ[1] * J₁)) * u_base(θ)
    du_2 = θ -> (exp(ψ[3] * J₃) * J₂ * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)
    du_3 = θ -> (J₃ * exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)

    du_11 = θ -> (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * J₁^2 * exp(ψ[1] * J₁)) * u_base(θ)
    du_12 = θ -> (exp(ψ[3] * J₃) * J₂ * exp(ψ[2] * J₂) * J₁ * exp(ψ[1] * J₁)) * u_base(θ)
    du_13 = θ -> (J₃ * exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * J₁ * exp(ψ[1] * J₁)) * u_base(θ)

    du_22 = θ -> (exp(ψ[3] * J₃) * J₂^2 * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)
    du_23 = θ -> (J₃ * exp(ψ[3] * J₃) * J₂ * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)

    du_33 = θ -> (J₃^2 * exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base(θ)

    # Computing derivatives with respect to ψ  

    derivative[1, 1] = c * fft_integrate(θ -> f_HA(u(θ), du_1(θ), du_1(θ), du_11(θ)) * interval(1) / (interval(1) + e * cos(θ)), interval(2) * π, N_fft, rigorous=true)
    derivative[1, 2] = c * fft_integrate(θ -> f_HA(u(θ), du_1(θ), du_2(θ), du_12(θ)) * 1interval(1) / (interval(1) + e * cos(θ)), interval(2) * π, N_fft, rigorous=true)
    derivative[1, 3] = c * fft_integrate(θ -> f_HA(u(θ), du_1(θ), du_3(θ), du_13(θ)) * interval(1) / (interval(1) + e * cos(θ)), interval(2) * π, N_fft, rigorous=true)

    derivative[2, 1] = derivative[1, 2]
    derivative[2, 2] = c * fft_integrate(θ -> f_HA(u(θ), du_2(θ), du_2(θ), du_22(θ)) * interval(1) / (interval(1) + e * cos(θ)), interval(2) * π, N_fft, rigorous=true)
    derivative[2, 3] = c * fft_integrate(θ -> f_HA(u(θ), du_2(θ), du_3(θ), du_23(θ)) * interval(1) / (interval(1) + e * cos(θ)), interval(2) * π, N_fft, rigorous=true)

    derivative[3, 1] = derivative[1, 3]
    derivative[3, 2] = derivative[2, 3]
    derivative[3, 3] = c * fft_integrate(θ -> f_HA(u(θ), du_3(θ), du_3(θ), du_33(θ)) * interval(1) / (interval(1) + e * cos(θ)), interval(2) * π, N_fft, rigorous=true)

    # Computing derivatives with respect to e 

    dc = c^4 / (interval(3) * pi) * fft_integrate(θ -> cos(θ) / ((interval(1) + e * cos(θ))^3), interval(2) * π, N_fft, rigorous=true)
    d2c = interval(4) * c^3 * dc / (interval(3) * pi) * fft_integrate(θ -> cos(θ) / ((interval(1) + e * cos(θ))^3), interval(2) * π, N_fft, rigorous=true) -
          c^4 / pi * fft_integrate(θ -> cos(θ)^2 / ((interval(1) + e * cos(θ))^4), interval(2) * π, N_fft, rigorous=true)

    de_f = θ -> dc / (interval(1) + e * cos(θ)) - c * cos(θ) / ((interval(1) + e * cos(θ))^2)
    de_2_f = θ -> d2c / (interval(1) + e * cos(θ)) -
                  interval(2) * dc * cos(θ) / ((interval(1) + e * cos(θ))^2) +
                  interval(2) * c * (cos(θ)^2) / ((interval(1) + e * cos(θ))^3)


    derivative[4, 4] = fft_integrate(θ -> f_A(u(θ)) * de_2_f(θ), interval(2) * π, N_fft, rigorous=true)

    # Computing mixed derivatives
    derivative[1, 4] = fft_integrate(θ -> f_DA(u(θ), du_1(θ)) * de_f(θ), interval(2) * π, N_fft, rigorous=true)
    derivative[4, 1] = derivative[1, 4]

    derivative[2, 4] = fft_integrate(θ -> f_DA(u(θ), du_2(θ)) * de_f(θ), interval(2) * π, N_fft, rigorous=true)
    derivative[4, 2] = derivative[2, 4]

    derivative[3, 4] = fft_integrate(θ -> f_DA(u(θ), du_3(θ)) * de_f(θ), interval(2) * π, N_fft, rigorous=true)
    derivative[4, 3] = derivative[3, 4]

    return derivative
end

# function f_A_vector(u, θ)
#     gens, num_gen = generators()
#     sum = zeros(length(θ), 1)

#     for i ∈ 1:num_gen
#         A = I - gens[:, :, i]
#         Au = A * u

#         sum = sum + 1 ./ ((Au[1, :] .^ 2 + Au[2, :] .^ 2 + Au[3, :] .^ 2) .^ (1 / 2))
#     end

#     return sum
# end

# function f_DA_vector(u, du, θ)
#     gens, num_gen = generators()
#     tmp = zeros(length(θ), 1)

#     for i ∈ 1:num_gen
#         A = I - gens[:, :, i]
#         Au = A * u
#         Adu = A * du

#         tmp = tmp - 1 ./ ((Au[1, :] .^ 2+Au[2, :] .^ 2+Au[3, :] .^ 2) .^ (3/2))[:] .* sum(Au .* Adu, dims=1)[:]
#     end

#     return tmp
# end

# function f_HA_vector(u, du_1, du_2, du_12, θ)
#     gens, num_gen = generators()
#     tmp = zeros(length(θ), 1)

#     for i ∈ 1:num_gen
#         A = I - gens[:, :, i]
#         Au = A * u
#         Adu_1 = A * du_1
#         Adu_2 = A * du_2
#         Adu_12 = A * du_12

#         tmp = tmp - 1 ./ ((Au[1, :] .^ 2+Au[2, :] .^ 2+Au[3, :] .^ 2) .^ (3/2))[:] .* sum(Adu_2 .* Adu_1, dims=1)[:]
#         tmp = tmp - 1 ./ ((Au[1, :] .^ 2+Au[2, :] .^ 2+Au[3, :] .^ 2) .^ (3/2))[:] .* sum(Au .* Adu_12, dims=1)[:]
#         tmp = tmp + 3 ./ ((Au[1, :] .^ 2+Au[2, :] .^ 2+Au[3, :] .^ 2) .^ (5/2))[:] .* sum(Au .* Adu_1, dims=1)[:] .* sum(Au .* Adu_2, dims=1)[:]
#     end

#     return tmp
# end

# function A_1_vector(X::Sequence, N_fft::Int64)
#     # Extracting input
#     ψ = component(X, 1)
#     e = component(X, 2)[1]

#     # Computing constant c
#     θ = collect(LinRange(0, 2 * pi, N_fft + 1))[1:end-1]
#     f_θ = (1 .+ e * cos.(θ)) .^ -2
#     c = (N_fft / FFTW.fft(f_θ)[1])^(1 / 3)

#     J₁ = [
#         0 0 0
#         0 0 -1
#         0 1 0
#     ]

#     J₂ = [
#         0 0 -1
#         0 0 0
#         1 0 0
#     ]

#     J₃ = [
#         0 -1 0
#         1 0 0
#         0 0 0
#     ]

#     rot = exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)

#     u = rot * [transpose(cos.(θ)); transpose(sin.(θ)); zeros(1, length(θ))]

#     θ = collect(LinRange(0, 2 * pi, N_fft + 1))[1:end-1]
#     f_θ = f_A_vector(u, θ) .* 1 ./ (1 .+ e * cos.(θ))

#     return c * 2 * pi / N_fft * FFTW.fft(f_θ)[1]
# end

# function DA_1_vector(X::Sequence, N_fft::Int64)
#     # Extracting input
#     ψ = component(X, 1)
#     e = component(X, 2)[1]

#     derivative = zeros(4, 1)
#     θ = collect(LinRange(0, 2 * pi, N_fft + 1))[1:end-1]

#     # Computing constant c
#     θ = collect(LinRange(0, 2 * pi, N_fft + 1))[1:end-1]
#     f_θ = (1 .+ e * cos.(θ)) .^ -2
#     c = (N_fft / FFTW.fft(f_θ)[1])^(1 / 3)

#     # Defining rotation matrices
#     J₁ = [
#         0 0 0
#         0 0 -1
#         0 1 0
#     ]

#     J₂ = [
#         0 0 -1
#         0 0 0
#         1 0 0
#     ]

#     J₃ = [
#         0 -1 0
#         1 0 0
#         0 0 0
#     ]

#     u_base = [transpose(cos.(θ)); transpose(sin.(θ)); zeros(1, length(θ))]

#     u = (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base

#     du_1 = (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * J₁ * exp(ψ[1] * J₁)) * u_base
#     du_2 = (exp(ψ[3] * J₃) * J₂ * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base
#     du_3 = (J₃ * exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base

#     # Computing derivatives with respect to ψ  

#     # ψ₁
#     derivative[1] = 2 * pi * c / N_fft * FFTW.fft(f_DA_vector(u, du_1, θ) .* 1 ./ (1 .+ e * cos.(θ)))[1]
#     # ψ₂
#     derivative[2] = 2 * pi * c / N_fft * FFTW.fft(f_DA_vector(u, du_2, θ) .* 1 ./ (1 .+ e * cos.(θ)))[1]
#     # ψ₃
#     derivative[3] = 2 * pi * c / N_fft * FFTW.fft(f_DA_vector(u, du_3, θ) .* 1 ./ (1 .+ e * cos.(θ)))[1]

#     # Computing derivative with respect to e 

#     dc = c^4 / (3 * pi) * 2 * pi / N_fft * FFTW.fft(cos.(θ) ./ ((1 .+ e * cos.(θ)) .^ 3))[1]
#     de_f = dc ./ (1 .+ e * cos.(θ)) - c * cos.(θ) ./ ((1 .+ e * cos.(θ)) .^ 2)
#     derivative[4] = 2 * pi / N_fft * FFTW.fft(f_A_vector(u, θ) .* de_f)[1]

#     return derivative
# end

# function HA_1_vector(X::Sequence, N_fft::Int64)
#     # Extracting input
#     ψ = X[1:3]
#     e = X[4]

#     derivative = zeros(4, 4)
#     θ = collect(LinRange(0, 2 * pi, N_fft + 1))[1:end-1]

#     # Computing constant c
#     θ = collect(LinRange(0, 2 * pi, N_fft + 1))[1:end-1]
#     f_θ = (1 .+ e * cos.(θ)) .^ -2
#     c = (N_fft / FFTW.fft(f_θ)[1])^(1 / 3)

#     # Defining rotation matrices
#     J₁ = [
#         0 0 0
#         0 0 -1
#         0 1 0
#     ]

#     J₂ = [
#         0 0 -1
#         0 0 0
#         1 0 0
#     ]

#     J₃ = [
#         0 -1 0
#         1 0 0
#         0 0 0
#     ]

#     u_base = [transpose(cos.(θ)); transpose(sin.(θ)); zeros(1, length(θ))]

#     u = (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base

#     du_1 = (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * J₁ * exp(ψ[1] * J₁)) * u_base
#     du_2 = (exp(ψ[3] * J₃) * J₂ * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base
#     du_3 = (J₃ * exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base

#     du_11 = (exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * J₁^2 * exp(ψ[1] * J₁)) * u_base
#     du_12 = (exp(ψ[3] * J₃) * J₂ * exp(ψ[2] * J₂) * J₁ * exp(ψ[1] * J₁)) * u_base
#     du_13 = (J₃ * exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * J₁ * exp(ψ[1] * J₁)) * u_base

#     du_22 = (exp(ψ[3] * J₃) * J₂^2 * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base
#     du_23 = (J₃ * exp(ψ[3] * J₃) * J₂ * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base

#     du_33 = (J₃^2 * exp(ψ[3] * J₃) * exp(ψ[2] * J₂) * exp(ψ[1] * J₁)) * u_base

#     # Computing derivatives with respect to ψ  

#     derivative[1, 1] = 2 * pi * c / N_fft * FFTW.fft(f_HA_vector(u, du_1, du_1, du_11, θ) .* 1 ./ (1 .+ e * cos.(θ)))[1]
#     derivative[1, 2] = 2 * pi * c / N_fft * FFTW.fft(f_HA_vector(u, du_1, du_2, du_12, θ) .* 1 ./ (1 .+ e * cos.(θ)))[1]
#     derivative[1, 3] = 2 * pi * c / N_fft * FFTW.fft(f_HA_vector(u, du_1, du_3, du_13, θ) .* 1 ./ (1 .+ e * cos.(θ)))[1]

#     derivative[2, 1] = derivative[1, 2]
#     derivative[2, 2] = 2 * pi * c / N_fft * FFTW.fft(f_HA_vector(u, du_2, du_2, du_22, θ) .* 1 ./ (1 .+ e * cos.(θ)))[1]
#     derivative[2, 3] = 2 * pi * c / N_fft * FFTW.fft(f_HA_vector(u, du_2, du_3, du_23, θ) .* 1 ./ (1 .+ e * cos.(θ)))[1]

#     derivative[3, 1] = derivative[1, 3]
#     derivative[3, 2] = derivative[2, 3]
#     derivative[3, 3] = 2 * pi * c / N_fft * FFTW.fft(f_HA_vector(u, du_3, du_3, du_33, θ) .* 1 ./ (1 .+ e * cos.(θ)))[1]

#     # Computing derivatives with respect to e 

#     dc = c^4 / (3 * pi) * 2 * pi / N_fft * FFTW.fft(cos.(θ) ./ ((1 .+ e * cos.(θ)) .^ 3))[1]
#     d2c = 4 * c^3 * dc / (3 * pi) * 2 * pi / N_fft * FFTW.fft(cos.(θ) ./ ((1 .+ e * cos.(θ)) .^ 3))[1] -
#           c^4 / pi * 2 * pi / N_fft * FFTW.fft((cos.(θ) .^ 2) ./ ((1 .+ e * cos.(θ)) .^ 4))[1]

#     de_f = dc ./ (1 .+ e * cos.(θ)) - c * cos.(θ) ./ ((1 .+ e * cos.(θ)) .^ 2)
#     de_2_f = d2c ./ (1 .+ e * cos.(θ)) -
#              2 * dc * cos.(θ) ./ ((1 .+ e * cos.(θ)) .^ 2) +
#              2 * c * (cos.(θ) .^ 2) ./ ((1 .+ e * cos.(θ)) .^ 3)


#     derivative[4, 4] = 2 * pi / N_fft * FFTW.fft(f_A_vector(u, θ) .* de_2_f)[1]

#     # Computing mixed derivatives
#     derivative[1, 4] = 2 * pi / N_fft * FFTW.fft(f_DA_vector(u, du_1, θ) .* de_f)[1]
#     derivative[4, 1] = derivative[1, 4]

#     derivative[2, 4] = 2 * pi / N_fft * FFTW.fft(f_DA_vector(u, du_2, θ) .* de_f)[1]
#     derivative[4, 2] = derivative[2, 4]

#     derivative[3, 4] = 2 * pi / N_fft * FFTW.fft(f_DA_vector(u, du_3, θ) .* de_f)[1]
#     derivative[4, 3] = derivative[3, 4]

#     return derivative
# end