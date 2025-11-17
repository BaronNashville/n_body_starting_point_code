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

# function A_1(ψ, e, N_fft)
#     sample_points, sample_time = kepler_sample(ψ, e, 0.0, N_fft)

#     function f_A(u)
#         gens, num_gen = generators()
#         sum = 0    

#         for i ∈ 1:num_gen
#             A = I - gens[:,:,i]
#             Au = A*u

#             # if norm(Au) < 1e-6
#             #     println("Particles are too close\nProblem occured using u =" * string(u) * "and generator g =" * string(gens[:,:,i]))
#             # end
#             sum = sum + 1/((Au[1]^2 + Au[2]^2 + Au[3]^2)^(1/2))
#             # println(sum)
#         end

#         return sum
#     end
#     f_points = zeros(1, N_fft)
#     for i ∈ 1:N_fft
#         f_points[i] = f_A(sample_points[:,i])
#     end

#     # display(f_points)

#     return 2*pi/N_fft * FFTW.fft(f_points)[1]
# end

function f_A(u, θ)
    gens, num_gen = generators()
    sum = zeros(length(θ),1)

    for i ∈ 1:num_gen
        A = I - gens[:,:,i]
        Au = A*u

        sum = sum + 1 ./((Au[1,:].^2 + Au[2,:].^2 + Au[3,:].^2).^(1/2))
    end

    return sum
end

function f_DA(u, du, θ)
    gens, num_gen = generators()
    tmp = zeros(length(θ),1)   

    for i ∈ 1:num_gen
        A = I - gens[:,:,i]
        Au = A*u
        Adu = A*du

        tmp = tmp - 1 ./((Au[1,:].^2 + Au[2,:].^2 + Au[3,:].^2).^(3/2)) .* sum(Au .* Adu, 1)
    end

    return tmp
end

function A_1(ψ, e, N_fft)
    # Compting constant c
    θ = collect(LinRange(0, 2*pi, N_fft+1))[1:end-1]
    f_θ = (1 .+e*cos.(θ)).^-2
    c = (N_fft/ FFTW.fft(f_θ)[1])^(1/3)

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

    u = rot * [transpose(cos.(θ)); transpose(sin.(θ)); zeros(1, length(θ))]

    θ = collect(LinRange(0,2*pi, N_fft+1))[1:end-1]
    f_θ = f_A(u, θ) .* 1 ./ (1 .+ e*cos.(θ))

    return c*2*pi/N_fft * FFTW.fft(f_θ)[1]
end

function DA_1(ψ, e, N_fft)
    derivative = zeros(1,3)
    θ = collect(LinSpace(0, 2*pi, N_fft+1))[1:end-1]

    # Compting constant c
    θ = collect(LinRange(0, 2*pi, N_fft+1))[1:end-1]
    f_θ = (1 .+e*cos.(θ)).^-2
    c = (N_fft/ FFTW.fft(f_θ)[1])^(1/3)

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

    u_base = [transpose(cos.(θ)), transpose(sin.(θ)), zeros(1, length(θ))]
    u = (exp(ψ[3]*J₃) * exp(ψ[2]*J₂) * exp(ψ[1]*J₁)) * u_base

    du_1 = (exp(ψ[3]*J₃) * exp(ψ[2]*J₂) * J₁ * exp(ψ[1]*J₁)) * u_base
    du_2 = (exp(ψ[3]*J₃) * J₂ * exp(ψ[2]*J₂) * exp(ψ[1]*J₁)) * u_base
    du_3 = (J₃ * exp(ψ[3]*J₃) * exp(ψ[2]*J₂) * exp(ψ[1]*J₁)) * u_base

    # Computing derivatives with respect to ψ  

    # ψ₁
    derivative[1] = 2*pi*c /N_fft * FFTW.fft(f_DA(u, du_1, θ) .* 1 ./ (1 .+ e*cos.(θ)))[1]
    # ψ₂
    derivative[2] = 2*pi*c /N_fft * FFTW.fft(f_DA(u, du_2, θ) .* 1 ./ (1 .+ e*cos.(θ)))[1] 
    # ψ₃
    derivative[3] = 2*pi*c /N_fft * FFTW.fft(f_DA(u, du_3, θ) .* 1 ./ (1 .+ e*cos.(θ)))[1]
    
    # Computing derivatives with respect to e 

    c_de = c^4/(3*pi) * 2*pi/N_fft * FFTW.fft(cos.(θ) ./ ((1 .+ e*cos.(θ)).^3))[1]
    de_f = c_de ./ (1 .+ e*cos.(θ)) - c * 1 ./((1 .+ e*cos.(θ)).^2)
    derivative[4] = 2*pi*c /N_fft * FFTW.fft(f_A(u, θ) .* de_f)[1]

    return derivative
end