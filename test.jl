include("action.jl")
include("algorithms.jl")
include("proofs.jl")
include("approx_derivatives.jl")

crit_points = [
    [
        1.7757222200167417
        0.02114435587723842
        0.20492589505849618
        1.3366368953365808e-7
    ],
    [
        0.2093202529695168
        2.936713216628474
        1.549200239859475
        1.0341282567446812e-7
    ],
    [
        1.7757222204642675
        3.162737111284318
        1.365870452903194
        1.3458978553254087e-7
    ],
    [
        1.775722220164424535370728410810336339763
        0.02114438947681569714850869983590603282479
        0.2049258880734805233140389228178926128987
        1.700237229263850965863518226152121868421e-08
    ],
    [
        4.94091780239843
        0.7583354177380563
        4.386255434244429
        0.8448356513624653
    ],
    [
        1.3422675047811548
        0.758335417738055
        1.8969298729351591
        0.84483565136246
    ],
    [
        1.3422675047811556
        5.524849889441536
        4.386255434244428
        0.8448356513624582
    ],
    [
        4.940917802398434
        5.524849889441528
        1.896929872935161
        0.8448356513624576
    ]
]

# setprecision(64)

N_fft = 2^15

# ψ = 2π * rand(3)
# e = 0

# X = rand(4)
# println(HA_1_float(X, N_fft, generators_12)-HA_1_approx(X, N_fft, generators_12))

# X = crit_points[1]
# X = [3.9269908169872414, -1.2789504890390249, 2.648622743286701, 0]


# eigen(HA_1_float(X, N_fft, generators_12))
X, _ = MyNewton(X -> (f_float(X[1:2], N_fft, generators_12), Df_float(X[1:2], N_fft, generators_12)), X[1:2], tol=1e-14, max_iter=10)

#println(norm(inv(Df(X[1:2], N_fft, generators_12))))

# println(Df(X[1:2], N_fft, generators_12) \ f(X[1:2], N_fft, generators_12))

display(first_order_radii_polymnomial(interval.(X[1:2]), x -> f_interval(x, N_fft, generators_12), x -> Df_interval(x, N_fft, generators_12), r_star = 1e-10))

# H = new_hessian(interval.(X[1:2]), N_fft, generators_12)

# DA_1_interval([interval.(X[1:2]); interval(0); interval(0)], N_fft, generators_12)

# display(DA_1_float([X[1:2];0.0;0.0], N_fft, generators_12)[1:2])
# display(HA_1_float([X[1:2];0.0;0.0], N_fft, generators_12)[1:2,1:2])
# display(inv(HA_1_float([X[1:2];0.0;0.0], N_fft, generators_12)[1:2,1:2]))
# display(inv(HA_1_float(X, N_fft, generators_12)))

# # Extracting input
# ψ = X[1:3]
# e = X[4]

# # Defining the correct integral
# integrate = f -> fft_integrate_float(f, 2, N_fft)

# derivative = zeros(4, 1)

# # Computing constant c
# # f_c = θ -> (1 + e * cos(θ))^(-2)
# # c = (2 * (π) / integrate(f_c))^(1 / 3)
# c = real(sqrt(complex(1-e^2)))
# dc = -e/c

# # Defining rotation matrices
# J₁, J₂, J₃ = rotation_matrices()

# 𝒥₁ = [
#     1 0 0
#     0 cos(ψ[1]) -sin(ψ[1])
#     0 sin(ψ[1]) cos(ψ[1])
# ]
# 𝒥₂ = [
#     cos(ψ[2]) 0 -sin(ψ[2])
#     0 1 0
#     sin(ψ[2]) 0 cos(ψ[2])
# ]
# 𝒥₃ = [
#     cos(ψ[3]) -sin(ψ[3]) 0
#     sin(ψ[3]) cos(ψ[3]) 0
#     0 0 1
# ]

# # 𝒥₁ = exp(J₁)^ψ[1]
# # 𝒥₂ = exp(J₂)^ψ[2]
# # 𝒥₃ = exp(J₃)^ψ[3]

# # if T <: Real
# #     𝒥₁ = real(𝒥₁)
# #     𝒥₂ = real(𝒥₂)
# #     𝒥₃ = real(𝒥₃)
# # end

# u_base = θ -> [cos(θ); sin(θ); 0]

# u = θ -> 𝒥₁ * 𝒥₂ * 𝒥₃ * u_base(θ)

# du_1 = θ -> J₁ * 𝒥₁  * 𝒥₂ * 𝒥₃ * u_base(θ)
# du_2 = θ -> 𝒥₁ * J₂ * 𝒥₂ * 𝒥₃ * u_base(θ)
# du_3 = θ -> 𝒥₁ * 𝒥₂ * J₃ * 𝒥₃ * u_base(θ)

# # Computing derivatives with respect to ψ  
# # ψ₁
# derivative[1] = c * integrate(θ -> f_DA_float(u(θ), du_1(θ), generators) * 1 / (1 + e * cos(θ)))
# # ψ₂
# derivative[2] = c * integrate(θ -> f_DA_float(u(θ), du_2(θ), generators) * 1 / (1 + e * cos(θ)))
# # ψ₃
# derivative[3] = c * integrate(θ -> f_DA_float(u(θ), du_3(θ), generators) * 1 / (1 + e * cos(θ)))

# # Computing derivative with respect to e 

# # dc = (1 / 3) * c^4 * integrate(θ -> cos(θ) / ((1 + e * cos(θ))^3)) / (π)
# # derivative[4] = integrate(θ -> (dc / (1 + e * cos(θ)) - c * cos(θ) / ((1 + e * cos(θ))^2)) * f_A_float(u(θ), generators))

# origin = rand(3)

# generators, num_gen = generators_60()

# particles = zeros(3, num_gen)

# for i ∈ 1:num_gen
#     particles[:,i] = generators[:,:,i] * origin
# end

# display(particles)

# for i ∈ 1:num_gen
#     for j ∈ i+1:num_gen
#         if norm(particles[:,i] - particles[:,j]) < 1e-5
#             println((i,j))
#         end
#     end
# end
# GLMakie.scatter(particles)

# # Newton approach
# X = [2*pi*rand(3);0.8]
# sol, success = MyNewton(X -> (DA_1_float(X, N_fft, generators_60), HA_1_float(X, N_fft, generators_60)), X, tol = 1e-12, max_iter = 20, verbose = true)

