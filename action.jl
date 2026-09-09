include("integrate.jl")
include("helpers.jl")

###################################################################
# Floating point implementation of the action and its derivatives #
###################################################################

function f_A_float(u::Vector{T}, gens, num_gen) where T
    tmp = zero(T)

    for i ∈ 2:num_gen
        L = gens[:, :, i]
        Au = u - L * u

        tmp = tmp + 1 / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^(1 / 2))
    end

    return tmp
end

function f_DA_float(u::Vector{T}, du::Vector{T}, gens, num_gen) where T
    tmp = zero(T)

    for i ∈ 2:num_gen
        L = gens[:, :, i]
        Au = u - L * u
        Adu = du - L * du

        tmp = tmp - 1 / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^(3 / 2)) * sum(Au .* Adu)
    end
    return tmp
end

function f_HA_float(u::Vector{T}, du_1::Vector{T}, du_2::Vector{T}, du_12::Vector{T}, gens, num_gen) where T
    tmp = zero(T)

    for i ∈ 2:num_gen
        L = gens[:, :, i]
        Au = u - L * u
        Adu_1 = du_1 - L * du_1
        Adu_2 = du_2 - L * du_2
        Adu_12 = du_12 - L * du_12

        tmp = tmp - 1 / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^(3 / 2)) * sum(Adu_2 .* Adu_1)
        tmp = tmp - 1 / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^(3 / 2)) * sum(Au .* Adu_12)
        tmp = tmp + 3 / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^(5 / 2)) * sum(Au .* Adu_1) * sum(Au .* Adu_2)
    end

    return tmp
end

function A_1_float(X::Vector{T}, N_fft, gens, num_gen) where T
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    # Defining the correct integral
    integrate = f -> fft_integrate_float(f, 2, N_fft)

    # Computing constant c
    c = real(sqrt(complex(1-e^2)))

    # Defining rotation matrices
    J₁, J₂, J₃ = rotation_generators()

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

    u_base = θ -> [cos(θ); sin(θ); 0]

    u = θ -> 𝒥₁ * 𝒥₂ * 𝒥₃ * u_base(θ)

    return c * integrate(θ -> f_A_float(u(θ), gens, num_gen) * 1 / (1 + e * cos(θ)))
end

function DA_1_float(X::Vector{T}, N_fft, gens, num_gen) where T
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    # Defining the correct integral
    integrate = f -> fft_integrate_float(f, 2, N_fft)

    derivative = zeros(T, 4, 1)

    # Computing constant c
    c = real(sqrt(complex(1-e^2)))
    dc = -e/c

    # Defining rotation matrices
    J₁, J₂, J₃ = rotation_generators()

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

    u_base = θ -> [cos(θ); sin(θ); 0]

    u = θ -> 𝒥₁ * 𝒥₂ * 𝒥₃ * u_base(θ)

    du_1 = θ -> J₁ * 𝒥₁  * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_2 = θ -> 𝒥₁ * J₂ * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_3 = θ -> 𝒥₁ * 𝒥₂ * J₃ * 𝒥₃ * u_base(θ)

    # Computing derivatives with respect to ψ  
    # ψ₁
    derivative[1] = c * integrate(θ -> f_DA_float(u(θ), du_1(θ), gens, num_gen) * 1 / (1 + e * cos(θ)))
    # ψ₂
    derivative[2] = c * integrate(θ -> f_DA_float(u(θ), du_2(θ), gens, num_gen) * 1 / (1 + e * cos(θ)))
    # ψ₃
    derivative[3] = c * integrate(θ -> f_DA_float(u(θ), du_3(θ), gens, num_gen) * 1 / (1 + e * cos(θ)))

    # Computing derivative with respect to e 

    # dc = (1 / 3) * c^4 * integrate(θ -> cos(θ) / ((1 + e * cos(θ))^3)) / (π)
    derivative[4] = integrate(θ -> (dc / (1 + e * cos(θ)) - c * cos(θ) / ((1 + e * cos(θ))^2)) * f_A_float(u(θ), gens, num_gen))

    return derivative
end

function HA_1_float(X::Vector{T}, N_fft, gens, num_gen) where T
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    # Defining the correct integral
    integrate = f -> fft_integrate_float(f, 2, N_fft)

    derivative = zeros(T, 4, 4)

    # Computing constant c
    c = real(sqrt(complex(1-e^2)))
    dc = -e/c
    d2c = -c^(-3)

    # Defining rotation matrices
    J₁, J₂, J₃ = rotation_generators()

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

    u_base = θ -> [cos(θ); sin(θ); 0]

    u = θ -> 𝒥₁ * 𝒥₂ * 𝒥₃ * u_base(θ)

    du_1 = θ -> J₁ * 𝒥₁  * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_2 = θ -> 𝒥₁ * J₂ * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_3 = θ -> 𝒥₁ * 𝒥₂ * J₃ * 𝒥₃ * u_base(θ)

    du_11 = θ -> J₁ * J₁ * 𝒥₁ * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_12 = θ -> J₁ * 𝒥₁ * J₂ * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_13 = θ -> J₁ * 𝒥₁ * 𝒥₂ * J₃ * 𝒥₃ * u_base(θ)

    du_22 = θ -> 𝒥₁ * J₂ * J₂ * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_23 = θ -> 𝒥₁ * J₂ * 𝒥₂ * J₃ * 𝒥₃* u_base(θ)

    du_33 = θ -> 𝒥₁ * 𝒥₂ *  J₃ * J₃ * 𝒥₃* u_base(θ)

    # Computing derivatives with respect to ψ  
    derivative[1, 1] = c * integrate(θ -> f_HA_float(u(θ), du_1(θ), du_1(θ), du_11(θ), gens, num_gen) * 1 / (1 + e * cos(θ)))
    derivative[1, 2] = c * integrate(θ -> f_HA_float(u(θ), du_1(θ), du_2(θ), du_12(θ), gens, num_gen) * 1 / (1 + e * cos(θ)))
    derivative[1, 3] = c * integrate(θ -> f_HA_float(u(θ), du_1(θ), du_3(θ), du_13(θ), gens, num_gen) * 1 / (1 + e * cos(θ)))

    derivative[2, 1] = derivative[1, 2]
    derivative[2, 2] = c * integrate(θ -> f_HA_float(u(θ), du_2(θ), du_2(θ), du_22(θ), gens, num_gen) * 1 / (1 + e * cos(θ)))
    derivative[2, 3] = c * integrate(θ -> f_HA_float(u(θ), du_2(θ), du_3(θ), du_23(θ), gens, num_gen) * 1 / (1 + e * cos(θ)))

    derivative[3, 1] = derivative[1, 3]
    derivative[3, 2] = derivative[2, 3]
    derivative[3, 3] = c * integrate(θ -> f_HA_float(u(θ), du_3(θ), du_3(θ), du_33(θ), gens, num_gen) * 1 / (1 + e * cos(θ)))

    # Computing derivatives with respect to e 

    # dc = (1 / 3) * c^4 * integrate(θ -> cos(θ) / ((1 + e * cos(θ))^3)) / (π)
    # d2c = c^3 * dc * (4 / 3) * integrate(θ -> cos(θ) / ((1 + e * cos(θ))^3)) / (π) - c^4 / (π) * integrate(θ -> cos(θ)^2 / ((1 + e * cos(θ))^4))

    de_f = θ -> dc / (1 + e * cos(θ)) - c * cos(θ) / ((1 + e * cos(θ))^2)
    de_2_f = θ -> d2c / (1 + e * cos(θ)) -
                  2 * dc * cos(θ) / ((1 + e * cos(θ))^2) +
                  2 * c * (cos(θ)^2) / ((1 + e * cos(θ))^3)


    derivative[4, 4] = integrate(θ -> f_A_float(u(θ), gens, num_gen) * de_2_f(θ))

    # Computing mixed derivatives
    derivative[1, 4] = integrate(θ -> f_DA_float(u(θ), du_1(θ), gens, num_gen) * de_f(θ))
    derivative[4, 1] = derivative[1, 4]

    derivative[2, 4] = integrate(θ -> f_DA_float(u(θ), du_2(θ), gens, num_gen) * de_f(θ))
    derivative[4, 2] = derivative[2, 4]

    derivative[3, 4] = integrate(θ -> f_DA_float(u(θ), du_3(θ), gens, num_gen) * de_f(θ))
    derivative[4, 3] = derivative[3, 4]

    return derivative
end

########################################################################
# Interval arithmetic implementation of the action and its derivatives #
########################################################################

function f_A_interval(u::Vector{T}, gens, num_gen) where T
    tmp = zero(T)

    for i ∈ 2:num_gen
        L = interval.(gens[:, :, i])
        Au = u - L * u

        tmp = tmp + interval(1) / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^(interval(1) / interval(2)))
    end

    return tmp
end

function f_DA_interval(u::Vector{T}, du::Vector{T}, gens, num_gen) where T
    tmp = zero(T)

    for i ∈ 2:num_gen
        L = interval.(gens[:, :, i])
        Au = u - L * u
        Adu = du - L * du

        tmp = tmp - interval(1) / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^(interval(3) / interval(2))) * sum(Au .* Adu)
    end

    return tmp
end

function f_HA_interval(u::Vector{T}, du_1::Vector{T}, du_2::Vector{T}, du_12::Vector{T}, gens, num_gen) where T
    tmp = zero(T)

    for i ∈ 2:num_gen
        L = interval.(gens[:, :, i])
        Au = u - L * u
        Adu_1 = du_1 - L * du_1
        Adu_2 = du_2 - L * du_2
        Adu_12 = du_12 - L * du_12

        tmp = tmp - interval(1) / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^(interval(3) / interval(2))) * sum(Adu_2 .* Adu_1)
        tmp = tmp - interval(1) / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^(interval(3) / interval(2))) * sum(Au .* Adu_12)
        tmp = tmp + interval(3) / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^(interval(5) / interval(2))) * sum(Au .* Adu_1) * sum(Au .* Adu_2)
    end

    return tmp
end

function A_1_interval(X::Vector{T}, N_fft::Int64, gens, num_gen) where T
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    # Defining the correct integral
    integrate = f -> fft_integrate(f, 2 // 1, N_fft, rigorous=true)

    # Computing constant c
    # f_c = θ -> (interval(1) + e * cos(θ))^interval(-2)
    # c = (interval(2) * (π / integrate(f_c)))^(interval(1) / interval(3))
    c = sqrt(interval(1) - e^2)

    # Defining rotation matrices
    J₁, J₂, J₃ = interval.(rotation_generators())

    i1 = interval(1)
    i0 = interval(0)

    𝒥₁ = interval.([
        i1 i0 i0
        i0 cos(ψ[1]) -sin(ψ[1])
        i0 sin(ψ[1]) cos(ψ[1])
    ])
    𝒥₂ = interval.([
        cos(ψ[2]) i0 -sin(ψ[2])
        i0 i1 i0
        sin(ψ[2]) i0 cos(ψ[2])
    ])
    𝒥₃ = interval.([
        cos(ψ[3]) -sin(ψ[3]) i0
        sin(ψ[3]) cos(ψ[3]) i0
        i0 i0 i1
    ])

    u_base = θ -> [cos(θ); sin(θ); i0]

    u = θ -> 𝒥₁ * 𝒥₂ * 𝒥₃ * u_base(θ)

    return c * integrate(θ -> f_A_interval(u(θ), gens, num_gen) * i1 / (i1 + e * cos(θ)))
end

function DA_1_interval(X::Vector{T}, N_fft::Int64, gens, num_gen) where T
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    # Defining the correct integral
    integrate = f -> fft_integrate(f, 2 // 1, N_fft, rigorous=true)

    derivative = zeros(T, 4, 1)

    # Computing constant c
    c = sqrt(interval(1) - e^2)
    dc = -e / c

    # Defining rotation matrices
    J₁, J₂, J₃ = interval.(rotation_generators())

    i1 = interval(1)
    i0 = interval(0)

    𝒥₁ = interval.([
        i1 i0 i0
        i0 cos(ψ[1]) -sin(ψ[1])
        i0 sin(ψ[1]) cos(ψ[1])
    ])
    𝒥₂ = interval.([
        cos(ψ[2]) i0 -sin(ψ[2])
        i0 i1 i0
        sin(ψ[2]) i0 cos(ψ[2])
    ])
    𝒥₃ = interval.([
        cos(ψ[3]) -sin(ψ[3]) i0
        sin(ψ[3]) cos(ψ[3]) i0
        i0 i0 i1
    ])

    u_base = θ -> [cos(θ); sin(θ); i0]

    u = θ -> 𝒥₁ * 𝒥₂ * 𝒥₃ * u_base(θ)

    du_1 = θ -> J₁ * 𝒥₁  * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_2 = θ -> 𝒥₁ * J₂ * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_3 = θ -> 𝒥₁ * 𝒥₂ * J₃ * 𝒥₃ * u_base(θ)

    # Computing derivatives with respect to ψ  
    # ψ₁
    # println("Computing derivative with respect to ψ₁")
    derivative[1] = c * integrate(θ -> f_DA_interval(u(θ), du_1(θ), gens, num_gen) * i1 / (i1 + e * cos(θ)))
    # ψ₂
    # println("Computing derivative with respect to ψ₂")
    derivative[2] = c * integrate(θ -> f_DA_interval(u(θ), du_2(θ), gens, num_gen) * i1 / (i1 + e * cos(θ)))
    # ψ₃
    # println("Computing derivative with respect to ψ₃")
    derivative[3] = c * integrate(θ -> f_DA_interval(u(θ), du_3(θ), gens, num_gen) * i1/ (i1 + e * cos(θ)))

    # Computing derivative with respect to e 
    # println("Computing derivative with respect to e")
    # dc = interval(1 // 3) * c^4 * integrate(θ -> cos(θ) / ((interval(1) + e * cos(θ))^3)) / π
    derivative[4] = integrate(θ -> (dc / (i1 + e * cos(θ)) - c * cos(θ) / ((i1 + e * cos(θ))^2)) * f_A_interval(u(θ), gens, num_gen))

    return derivative
end

function HA_1_interval(X::Vector{T}, N_fft::Int64, gens, num_gen) where T
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    # Defining the correct integral
    integrate = f -> fft_integrate(f, 2 // 1, N_fft, rigorous=true)

    derivative = zeros(T, 4, 4)

    # Computing constant c
    c = sqrt(interval(1) - e^2)
    dc = -e/c
    d2c = -c^interval(-3)

    # Defining rotation matrices
    J₁, J₂, J₃ = interval.(rotation_generators())

    i1 = interval(1)
    i0 = interval(0)

    𝒥₁ = interval.([
        i1 i0 i0
        i0 cos(ψ[1]) -sin(ψ[1])
        i0 sin(ψ[1]) cos(ψ[1])
    ])
    𝒥₂ = interval.([
        cos(ψ[2]) i0 -sin(ψ[2])
        i0 i1 i0
        sin(ψ[2]) i0 cos(ψ[2])
    ])
    𝒥₃ = interval.([
        cos(ψ[3]) -sin(ψ[3]) i0
        sin(ψ[3]) cos(ψ[3]) i0
        i0 i0 i1
    ])

    u_base = θ -> [cos(θ); sin(θ); i0]

    u = θ -> 𝒥₁ * 𝒥₂ * 𝒥₃ * u_base(θ)

    du_1 = θ -> J₁ * 𝒥₁  * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_2 = θ -> 𝒥₁ * J₂ * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_3 = θ -> 𝒥₁ * 𝒥₂ * J₃ * 𝒥₃ * u_base(θ)

    du_11 = θ -> J₁ * J₁ * 𝒥₁ * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_12 = θ -> J₁ * 𝒥₁ * J₂ * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_13 = θ -> J₁ * 𝒥₁ * 𝒥₂ * J₃ * 𝒥₃ * u_base(θ)

    du_22 = θ -> 𝒥₁ * J₂ * J₂ * 𝒥₂ * 𝒥₃ * u_base(θ)
    du_23 = θ -> 𝒥₁ * J₂ * 𝒥₂ * J₃ * 𝒥₃* u_base(θ)

    du_33 = θ -> 𝒥₁ * 𝒥₂ *  J₃ * J₃ * 𝒥₃* u_base(θ)

    # Computing derivatives with respect to ψ  

    derivative[1, 1] = c * integrate(θ -> f_HA_interval(u(θ), du_1(θ), du_1(θ), du_11(θ), gens, num_gen) * i1 / (i1 + e * cos(θ)))
    derivative[1, 2] = c * integrate(θ -> f_HA_interval(u(θ), du_1(θ), du_2(θ), du_12(θ), gens, num_gen) * i1 / (i1 + e * cos(θ)))
    derivative[1, 3] = c * integrate(θ -> f_HA_interval(u(θ), du_1(θ), du_3(θ), du_13(θ), gens, num_gen) * i1 / (i1 + e * cos(θ)))

    derivative[2, 1] = derivative[1, 2]
    derivative[2, 2] = c * integrate(θ -> f_HA_interval(u(θ), du_2(θ), du_2(θ), du_22(θ), gens, num_gen) * i1 / (i1 + e * cos(θ)))
    derivative[2, 3] = c * integrate(θ -> f_HA_interval(u(θ), du_2(θ), du_3(θ), du_23(θ), gens, num_gen) * i1 / (i1 + e * cos(θ)))

    derivative[3, 1] = derivative[1, 3]
    derivative[3, 2] = derivative[2, 3]
    derivative[3, 3] = c * integrate(θ -> f_HA_interval(u(θ), du_3(θ), du_3(θ), du_33(θ), gens, num_gen) * i1 / (i1 + e * cos(θ)))

    # Computing derivatives with respect to e 
    de_f = θ -> dc / (i1 + e * cos(θ)) - c * cos(θ) / ((i1 + e * cos(θ))^2)
    de_2_f = θ -> d2c / (i1 + e * cos(θ)) -
                  interval(2) * dc * cos(θ) / ((i1 + e * cos(θ))^2) +
                  interval(2) * c * (cos(θ)^2) / ((i1 + e * cos(θ))^3)


    derivative[4, 4] = integrate(θ -> f_A_interval(u(θ), gens, num_gen) * de_2_f(θ))

    # Computing mixed derivatives
    derivative[1, 4] = integrate(θ -> f_DA_interval(u(θ), du_1(θ), gens, num_gen) * de_f(θ))
    derivative[4, 1] = derivative[1, 4]

    derivative[2, 4] = integrate(θ -> f_DA_interval(u(θ), du_2(θ), gens, num_gen) * de_f(θ))
    derivative[4, 2] = derivative[2, 4]

    derivative[3, 4] = integrate(θ -> f_DA_interval(u(θ), du_3(θ), gens, num_gen) * de_f(θ))
    derivative[4, 3] = derivative[3, 4]

    return derivative
end

######################
# Special case e = 0 #
######################

function f_float(X::Vector{T}, N_fft::Int64, gens, num_gen) where T
    return DA_1_float([X;0;0], N_fft, gens, num_gen)[1:2]
end

function Df_float(X::Vector{T}, N_fft::Int64, gens, num_gen) where T
    return HA_1_float([X;0;0], N_fft, gens, num_gen)[1:2,1:2]
end

function f_interval(X::Vector{T}, N_fft::Int64, gens, num_gen) where T
    return DA_1_interval([X;interval(0);interval(0)], N_fft, gens, num_gen)[1:2]
end

function Df_interval(X::Vector{T}, N_fft::Int64, gens, num_gen) where T
    return HA_1_interval([X;interval(0);interval(0)], N_fft, gens, num_gen)[1:2,1:2]
end

function new_hessian(X::Vector{T}, N_fft::Int64, gens, num_gen) where T
    i0 = interval(0)
    i1 = interval(1)
    i2 = interval(2)
    i4 = interval(4)
    iπ = interval(π)

    derivative = zeros(T, 4, 4)

    H00 = HA_1_interval([X;i0;i0], N_fft, gens, num_gen)
    Hπ2 = HA_1_interval([X;iπ/i2;i0], N_fft, gens, num_gen)
    Hπ4 = HA_1_interval([X;iπ/i4;i0], N_fft, gens, num_gen)

    A = H00[4,4]
    B = Hπ2[4,4]
    C = Hπ4[4,4] - i1/i2*(A+B)

    # Normal derivatives
    derivative[1:2,1:2] = H00[1:2,1:2]
    
    derivative[3:4,3:4] = [
        A C
        C B
    ]

    # Mixed derivatives
    derivative[1,3] = H00[1,4]
    derivative[1,4] = Hπ2[1,4]
    derivative[2,3] = H00[2,4]
    derivative[2,4] = Hπ2[2,4]

    derivative[3:4,1:2] = transpose(derivative[1:2,3:4])

    return derivative
end