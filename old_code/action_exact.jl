include("integrate.jl")
include("helpers.jl")

function f_A(u::Vector{T}) where T
    gens, num_gen = generators()
    tmp = zero(T)

    for i ∈ 1:num_gen
        L = gens[:, :, i]
        Au = u - L * u

        tmp = tmp + exact(1) / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^exact(1 // 2))
    end

    return tmp
end

function f_DA(u::Vector{T}, du::Vector{T}) where T
    gens, num_gen = generators()
    tmp = zero(T)

    for i ∈ 1:num_gen
        L = gens[:, :, i]
        Au = u - L * u
        Adu = du - L * du

        tmp = tmp - exact(1) / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^exact(3 // 2)) * sum(Au .* Adu)
    end
    return tmp
end

function f_HA(u::Vector{T}, du_1::Vector{T}, du_2::Vector{T}, du_12::Vector{T}) where T
    gens, num_gen = generators()
    tmp = zero(T)

    for i ∈ 1:num_gen
        L = gens[:, :, i]
        Au = u - L * u
        Adu_1 = du_1 - L * du_1
        Adu_2 = du_2 - L * du_2
        Adu_12 = du_12 - L * du_12

        tmp = tmp - exact(1) / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^exact(3 // 2)) * sum(Adu_2 .* Adu_1)
        tmp = tmp - exact(1) / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^exact(3 // 2)) * sum(Au .* Adu_12)
        tmp = tmp + exact(3) / ((Au[1] * Au[1] + Au[2] * Au[2] + Au[3] * Au[3])^exact(5 // 2)) * sum(Au .* Adu_1) * sum(Au .* Adu_2)
    end

    return tmp
end


function A_1(X::Vector{T}, N_fft::Int64) where T
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    # Defining the correct integral
    if T <: Interval || T <: Complex{Interval{V}} where V<:Real
        integrate = f -> fft_integrate(f, 2 // 1, N_fft, rigorous=true)
    else
        integrate = f -> fft_integrate(f, 2 // 1, N_fft)
    end

    # Computing constant c
    f_c = θ -> (exact(1) + e * cos(θ))^exact(-2)
    c = (exact(2) * (π / integrate(f_c)))^exact(1 // 3)

    # Defining rotation matrices
    J₁, J₂, J₃ = rotation_matrices()

    𝒥₁ = exp(ψ[1] * J₁)
    𝒥₂ = exp(ψ[2] * J₂)
    𝒥₃ = exp(ψ[3] * J₃)

    if T <: Real
        𝒥₁ = real(𝒥₁)
        𝒥₂ = real(𝒥₂)
        𝒥₃ = real(𝒥₃)
    end

    u_base = θ -> [cos(θ); sin(θ); exact(0)]
    u = θ -> 𝒥₃ * 𝒥₂ * 𝒥₁ * u_base(θ)

    return c * integrate(θ -> f_A(u(θ)) * exact(1) / (exact(1) + e * cos(θ)))
end

function DA_1(X::Vector{T}, N_fft::Int64) where T
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    # Defining the correct integral
    if T <: Interval || T <: Complex{Interval{V}} where V<:Real
        integrate = f -> fft_integrate(f, 2 // 1, N_fft, rigorous=true)
    else
        integrate = f -> fft_integrate(f, 2 // 1, N_fft)
    end

    derivative = zeros(T, 4, 1)

    # Computing constant c
    f_c = θ -> (exact(1) + e * cos(θ))^exact(-2)
    c = (exact(2) * (π / integrate(f_c)))^exact(1 // 3)

    # Defining rotation matrices
    J₁, J₂, J₃ = rotation_matrices()

    𝒥₁ = exp(ψ[1] * J₁)
    𝒥₂ = exp(ψ[2] * J₂)
    𝒥₃ = exp(ψ[3] * J₃)

    if T <: Real
        𝒥₁ = real(𝒥₁)
        𝒥₂ = real(𝒥₂)
        𝒥₃ = real(𝒥₃)
    end

    u_base = θ -> [cos(θ); sin(θ); exact(0)]

    u = θ -> 𝒥₃ * 𝒥₂ * 𝒥₁ * u_base(θ)

    du_1 = θ -> 𝒥₃ * 𝒥₂ * J₁ * 𝒥₁ * u_base(θ)
    du_2 = θ -> 𝒥₃ * J₂ * 𝒥₂ * 𝒥₁ * u_base(θ)
    du_3 = θ -> J₃ * 𝒥₃ * 𝒥₂ * 𝒥₁ * u_base(θ)

    # Computing derivatives with respect to ψ  
    # ψ₁
    derivative[1] = c * integrate(θ -> f_DA(u(θ), du_1(θ)) * exact(1) / (exact(1) + e * cos(θ)))
    # ψ₂
    derivative[2] = c * integrate(θ -> f_DA(u(θ), du_2(θ)) * exact(1) / (exact(1) + e * cos(θ)))
    # ψ₃
    derivative[3] = c * integrate(θ -> f_DA(u(θ), du_3(θ)) * exact(1) / (exact(1) + e * cos(θ)))

    # Computing derivative with respect to e 

    dc = exact(1 // 3) * c^4 * integrate(θ -> cos(θ) / ((exact(1) + e * cos(θ))^3)) / π
    derivative[4] = integrate(θ -> (dc / (exact(1) + e * cos(θ)) - c * cos(θ) / ((exact(1) + e * cos(θ))^2)) * f_A(u(θ)))

    return derivative
end

function HA_1(X::Vector{T}, N_fft::Int64) where T
    # Extracting input
    ψ = X[1:3]
    e = X[4]

    # Defining the correct integral
    if T <: Interval || T <: Complex{Interval{V}} where V<:Real
        integrate = f -> fft_integrate(f, 2 // 1, N_fft, rigorous=true)
    else
        integrate = f -> fft_integrate(f, 2 // 1, N_fft)
    end

    derivative = zeros(T, 4, 4)

    # Computing constant c
    f_c = θ -> (exact(1) + e * cos(θ))^exact(-2)
    c = (exact(2) * (π / integrate(f_c)))^exact(1 // 3)

    # Defining rotation matrices
    J₁, J₂, J₃ = rotation_matrices()

    𝒥₁ = exp(ψ[1] * J₁)
    𝒥₂ = exp(ψ[2] * J₂)
    𝒥₃ = exp(ψ[3] * J₃)

    if T <: Real
        𝒥₁ = real(𝒥₁)
        𝒥₂ = real(𝒥₂)
        𝒥₃ = real(𝒥₃)
    end

    u_base = θ -> [cos(θ); sin(θ); exact(0)]

    u = θ -> 𝒥₃ * 𝒥₂ * 𝒥₁ * u_base(θ)

    du_1 = θ -> 𝒥₃ * 𝒥₂ * J₁ * 𝒥₁ * u_base(θ)
    du_2 = θ -> 𝒥₃ * J₂ * 𝒥₂ * 𝒥₁ * u_base(θ)
    du_3 = θ -> J₃ * 𝒥₃ * 𝒥₂ * 𝒥₁ * u_base(θ)

    du_11 = θ -> 𝒥₃ * 𝒥₂ * J₁ * (J₁ * 𝒥₁) * u_base(θ)
    du_12 = θ -> 𝒥₃ * J₂ * 𝒥₂ * J₁ * 𝒥₁ * u_base(θ)
    du_13 = θ -> J₃ * 𝒥₃ * 𝒥₂ * J₁ * 𝒥₁ * u_base(θ)

    du_22 = θ -> 𝒥₃ * J₂ * (J₂ * 𝒥₂) * 𝒥₁ * u_base(θ)
    du_23 = θ -> J₃ * 𝒥₃ * J₂ * 𝒥₂ * 𝒥₁ * u_base(θ)

    du_33 = θ -> J₃ * (J₃ * 𝒥₃) * 𝒥₂ * 𝒥₁ * u_base(θ)

    # Computing derivatives with respect to ψ  

    derivative[1, 1] = c * integrate(θ -> f_HA(u(θ), du_1(θ), du_1(θ), du_11(θ)) * exact(1) / (exact(1) + e * cos(θ)))
    derivative[1, 2] = c * integrate(θ -> f_HA(u(θ), du_1(θ), du_2(θ), du_12(θ)) * exact(1) / (exact(1) + e * cos(θ)))
    derivative[1, 3] = c * integrate(θ -> f_HA(u(θ), du_1(θ), du_3(θ), du_13(θ)) * exact(1) / (exact(1) + e * cos(θ)))

    derivative[2, 1] = derivative[1, 2]
    derivative[2, 2] = c * integrate(θ -> f_HA(u(θ), du_2(θ), du_2(θ), du_22(θ)) * exact(1) / (exact(1) + e * cos(θ)))
    derivative[2, 3] = c * integrate(θ -> f_HA(u(θ), du_2(θ), du_3(θ), du_23(θ)) * exact(1) / (exact(1) + e * cos(θ)))

    derivative[3, 1] = derivative[1, 3]
    derivative[3, 2] = derivative[2, 3]
    derivative[3, 3] = c * integrate(θ -> f_HA(u(θ), du_3(θ), du_3(θ), du_33(θ)) * exact(1) / (exact(1) + e * cos(θ)))

    # Computing derivatives with respect to e 

    dc = exact(1 // 3) * c^4 * integrate(θ -> cos(θ) / ((exact(1) + e * cos(θ))^3)) / π
    d2c = c^3 * dc * exact(4 // 3) * integrate(θ -> cos(θ) / ((exact(1) + e * cos(θ))^3)) / π - c^4 / pi * integrate(θ -> cos(θ)^2 / ((exact(1) + e * cos(θ))^4))

    de_f = θ -> dc / (exact(1) + e * cos(θ)) - c * cos(θ) / ((exact(1) + e * cos(θ))^2)
    de_2_f = θ -> d2c / (exact(1) + e * cos(θ)) -
                  exact(2) * dc * cos(θ) / ((exact(1) + e * cos(θ))^2) +
                  exact(2) * c * (cos(θ)^2) / ((exact(1) + e * cos(θ))^3)


    derivative[4, 4] = integrate(θ -> f_A(u(θ)) * de_2_f(θ))

    # Computing mixed derivatives
    derivative[1, 4] = integrate(θ -> f_DA(u(θ), du_1(θ)) * de_f(θ))
    derivative[4, 1] = derivative[1, 4]

    derivative[2, 4] = integrate(θ -> f_DA(u(θ), du_2(θ)) * de_f(θ))
    derivative[4, 2] = derivative[2, 4]

    derivative[3, 4] = integrate(θ -> f_DA(u(θ), du_3(θ)) * de_f(θ))
    derivative[4, 3] = derivative[3, 4]

    return derivative
end