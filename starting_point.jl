include("action.jl")
include("proofs.jl")
include("algorithms.jl")

# We verify the existence of a non degenerate critical point

# crit_points_12 = [
#     [
#         2.35619   
#         1.33605
#         3.14159
#          0.844836
#     ],
#     [
#         2.35619
#         4.94713
#         3.14159
#         0.844836
#     ],
#     [
#         3.92699
#         4.94713
#         3.14159
#         0.844836
#     ]
# ]

# crit_points_60 = [
#     [
#         1.51207 
#         1.5708
#         4.72303
#         0.909775
#     ],
#     [
#         4.65366
#         1.5708
#         1.58144
#         0.909775
#     ],
#     [
#         4.37119
#         0.85582
#         2.94803
#         0.174288
#     ],
#     [
#         4.37119
#         2.28577
#         6.08962
#         -0.174288
#     ],
#     [
#         7.91271
#         4.71239
#         4.72303
#         0.909775
#     ],
#     [
#         5.05359
#         3.99741
#         2.94803
#         0.174288
#     ]
# ]

N_fft = 2^14

#########################
# Non-zero eccentricity #
#########################

# generators = generators_12
# X = [5.497787144; 4.947130974; 3.141592654; 0.844835651500]

# generators = generators_24
# X = [5.348741076; 5.200484038; 4.10039951; 0.873524144336]
# X = [5.441475127; 5.535368838; 2.949509909; 0.962477021261]

# generators = generators_60
# X = [5.793233620; 6.190758148; 6.259756780; 0.679505700484]

# X, success = MyNewton(X -> (DA_1_float(X, N_fft, generators), HA_1_float(X, N_fft, generators)), X, tol = 5e-14, max_iter = 4, verbose = true)

# (r, _, success) = first_order_radii_polymnomial(X, x -> DA_1_interval(x, N_fft, generators), x -> HA_1_interval(x, N_fft, generators), r_star = 1e-10)

# display((r, success))

#####################
# Zero eccentricity #
#####################

# generators = generators_12
# X = [5.497787144; 5.004234818]
# X = [5.497787144; 5.667705599]

# generators = generators_24
# X = [5.348522624; 5.503985790]
# X = [5.250347455; 5.166846112]

# generators = generators_60
# X = [5.799016554; 6.188895248]
# X = [5.849763501; 6.046873639]
# X = [5.983287023; 6.046539075]
# X = [6.191859601; 5.915899290]

X, success = MyNewton(X -> (f_float(X, N_fft, generators), Df_float(X, N_fft, generators)), X, tol = 5e-14, max_iter = 4, verbose = true)

(r, _, success) = first_order_radii_polymnomial(X, x -> f_interval(x, N_fft, generators), x -> Df_interval(x, N_fft, generators), r_star = 1e-10)

H = new_hessian(interval.(X).+interval(r)*interval(-1,1), N_fft, generators)

println("\ndet(H) = ")
display(det(H))


# nframes = 420
# framerate = 60
# time_step = 1:nframes

# X, _ = MyNewton(X -> (DA_1_float(X, N_fft, generators_60), HA_1_float(X, N_fft, generators_60)), crit_points_60[1])

# space_data, time_data = kepler_sample(X[1:3], X[4], 0.0, nframes)

# gens, num_gen = generators_60()

# points = Observable([Point3f(0,0,0) for _ in 1:num_gen])
# trails = [Observable{Vector{Point3f}}(Point3f[]) for _ in 1:num_gen]
# # colors = [:black, :red, :blue, :green, :orange,
# #           :purple, :cyan, :magenta, :yellow,
# #            :pink, :brown, :gray]


# fig = Figure(size = (1600,900))
# ax = Axis3(fig[1,1], title = L"Critical point of $A_1$", 
# 	  titlesize = 20, 
# 	  xlabel = L"$x$",
# 	  xlabelsize = 20,
# 	  ylabel = L"$y$",
# 	  ylabelsize = 20,
# 	  zlabel = L"$z$",
# 	  zlabelsize = 20,
#       limits=((-1.5,1.5), (-1.5,1.5), (-1.5,1.5))
# )

# scatter!(ax, points, #color = colors
# )
# scatter!(ax, Point3f(0,0,0), marker = :rect, color = :black, label = "Origin")

# for i in 1:num_gen
#     lines!(ax, trails[i], linewidth = 2, #color = (colors[i], 0.4)
#     )
# end

# axislegend("Legend")

# function change_function(time_step)
#     gens, num_gen = generators_60()
#     new_points = [Point3f(0,0,0) for _ in 1:num_gen]
#     for i ∈ 1:num_gen
#         p = Point3f(gens[:,:,i] * space_data[:,time_step])

#         new_points[i] = p

#         push!(trails[i][], p)
#         notify(trails[i])
#     end
#     # p0 = Point3f(space_data[:,time_step])

#     # new_points[12] = p0

#     # push!(trails[num_gen+1][], p0)
#     # notify(trails[num_gen+1])

#     points[] = new_points
# end

# record(change_function, fig, "crit_point_test.mp4", time_step; framerate = framerate)