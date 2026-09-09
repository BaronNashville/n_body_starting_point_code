include("action.jl")
include("proofs.jl")
include("algorithms.jl")

# We verify the existence of a non degenerate critical point of A_1
function main()
    ##############
    # PARAMETERS #
    ##############
    N_fft = 2^14

    proof = false

    movie = false
    picture = false

    __display__ = false
    __save__ = false
    __save_location__ = "./figures"

    show_duals = false
    num_duals = 1

    primal_outline_opacity = 1
    primal_mesh_opacity = 1
    primal_axes_opacity = 0.2

    if show_duals
        dual_outline_opacity = 0.1
        dual_mesh_opacity = 0.1
        dual_axes_opacity = 0
    else
        dual_outline_opacity = 0
        dual_mesh_opacity = 0
        dual_axes_opacity = 0
    end

    padding = 1.1

    ########
    # DATA #
    ########

    picture_color_pool = distinguishable_colors(22);
    movie_color_pool = [:black; :black; :red; fill(:blue, 19)];   

    #############
    # Eccentric #
    #############

    # DATA SHAPE : [ψ₁; ψ₂; ψ₃; e], (generators, number of generators, symmetries), name of saved files

    T_ecc_1 = [5.49778714378214; 4.947130974207735; 3.1415926535897944; 0.8448356513624544], generators_12(symmetries = true), "T_ecc_1"

    O_ecc_1 = [5.348741076103105; 5.200484037773405; 4.10039951058651; 0.8735241442595474], generators_24(symmetries = true), "O_ecc_1"
    O_ecc_2 = [5.441475127816426; 5.535368838122984; 2.949509909029401; 0.9624770211566274], generators_24(symmetries = true), "O_ecc_2"

    I_ecc_1 = [5.793233619038003; 6.190758148847939; 6.25975674713244; 0.6795057169690465], generators_60(symmetries = true), "I_ecc_1"

    ############
    # Circular #
    ############

    # DATA SHAPE : [ϵ₁; ϵ₂] , (generators, number of generators, symmetries), name of saved files

    T_circ_1 = [5.497787143782138; 5.004234818140562], generators_12(symmetries = true),"T_circ_1" # 8.0e2
    T_circ_2 = [5.497787143782138; 5.667705598509199], generators_12(symmetries = true), "T_circ_2" # none

    O_circ_1 = [5.348522623126767; 5.503985790434791], generators_24(symmetries = true), "O_circ_1" # 2.8e5
    O_circ_2 = [5.250347454029935; 5.166846111999057], generators_24(symmetries = true), "O_circ_2" # 1.0e5

    I_circ_1 = [5.7990165538259895; 6.188895248386591], generators_60(symmetries = true), "I_circ_1" # 1.15e7
    I_circ_2 = [5.84976350064333; 6.04687363867634], generators_60(symmetries = true), "I_circ_2" # 4.23e7
    I_circ_3 = [5.9832870230975175; 6.046539074513424], generators_60(symmetries = true), "I_circ_3" # 1.51e8
    I_circ_4 = [6.191859600880807; 5.915899290468376], generators_60(symmetries = true),"I_circ_4" # 1.1e8

    ################
    # All together #
    ################

    DATA = [T_ecc_1, O_ecc_1, O_ecc_2, I_ecc_1, T_circ_1, O_circ_1, O_circ_2, I_circ_1, I_circ_2, I_circ_3, I_circ_4]
    TEST = [O_ecc_1]

    for d ∈ DATA

        X, (gens, num_gen, syms), name = d

        eccentric = (length(X) == 4)

        if proof
            if eccentric
                (Y, Z, r, _, success) = first_order_radii_polymnomial(X, x -> DA_1_interval(x, N_fft, gens, num_gen), x -> HA_1_interval(x, N_fft, gens, num_gen), r_star = 1e-10, verbose = false)

                if success
                    println("\n" * name * " proof successful\nY: " * string(Y), "\nZ: " * string(Z) * "\nradius: " * string(r))
                else
                    println("\n" * name * " radii polynomial failed")
                end
            else
                (Y, Z, r, _, success) = first_order_radii_polymnomial(X, x -> f_interval(x, N_fft, gens, num_gen), x -> Df_interval(x, N_fft, gens, num_gen), r_star = 1e-10, verbose = false)


                H = new_hessian(interval.(X).+interval(r)*interval(-1,1), N_fft, gens, num_gen)
                detH = inf(abs(det(H)))

                if success && detH > 0
                    println("\n" * name * " proof successful\nY: " * string(Y), "\nZ: " * string(Z) * "\nradius: " * string(r) * "\nHessian determinant ≥ " * string(detH))
                elseif !success
                    println("\n" * name * " radii polynomial failed")
                elseif detH <= 0
                    println("\n" * name * " Hessian might be non invertible")
                end
            end            
        end
        
        if eccentric
            ψ = X[1:3]
            e = X[4]
        else
            ψ = [X[1:2];0]
            e = 0.0
        end

        generating_orbit = kepler_shape(ψ, e, 0.0, 1000)
        M = maximum(sum(generating_orbit.^2, dims=1).^(1/2))
        rot = rotation_matrix(ψ)

        if movie

            nframes = 1000
            framerate = 60
            time_step = 1:nframes

            s_data = kepler_shape(ψ, e, 3.14159, nframes)

            points = Observable([Point3f(0,0,0) for _ in 1:num_gen])

            fig = Figure(size = (800,800), padding = 0)
            ax = Axis3(fig[1,1],
                perspectiveness = 0.5,
                azimuth = 1.3 * π,
                elevation = π/8,
                aspect = :data,
                titlesize = 20, 
                xlabel = L"$x$",
                xlabelsize = 20,
                ylabel = L"$y$",
                ylabelsize = 20,
                zlabel = L"$z$",
                zlabelsize = 20,
                limits = ((-padding*M,padding*M), (-padding*M,padding*M), (-padding*M,padding*M))
            )

            outline_plots = []
            axes_plots = []
            mesh_plots = []

            syms_color = []

            for sym ∈ syms
                outline_sym_plots = []
                axes_sym_plots = []
                mesh_sym_plots = []

                push!(outline_plots, outline_sym_plots)
                push!(axes_plots, axes_sym_plots)
                push!(mesh_plots, mesh_sym_plots)

                # Extracting symmetry structure
                shape_sides, num_shapes = size(sym)
                sym_colors = repeat(movie_color_pool[3:num_shapes+2], inner=shape_sides)
                push!(syms_color, sym_colors)

                # Plotting the mesh and outline
                for i = 1:num_shapes
                    starting_index = 1 + (i-1)*shape_sides

                    polygon_points = @lift begin
                        p = Vector{Point3f}(undef, shape_sides)
                        for j in 1:shape_sides
                            p[j] = $(points)[sym[starting_index + j - 1]]
                        end
                        p
                    end

                    polygon_outline = @lift begin
                        p = Vector{Point3f}(undef, shape_sides + 1)
                        for j in 1:shape_sides
                            p[j] = $(points)[sym[starting_index + j - 1]]
                        end
                        p[end] = $(points)[sym[starting_index]]
                        p
                    end

                    # Outline            
                    push!(outline_sym_plots, lines!(ax, polygon_outline, color = sym_colors[starting_index]))

                    # Filled face
                    push!(mesh_sym_plots, mesh!(ax, polygon_points, triangulation(shape_sides), color = (sym_colors[starting_index
                    ], dual_mesh_opacity)))
                end

                # Plotting the rotation axes
                for i = 1:num_shapes
                    starting_index = 1 + (i-1)*shape_sides
                    vertices = zeros(3, shape_sides)
                    for j = 1:shape_sides
                        current_index = starting_index + j - 1
                        vertices[:,j] = gens[:,:,sym[current_index]] * s_data[:,Int(nframes/2)]
                    end
                    average = 1 / shape_sides * sum(vertices, dims = 2)
                    rot_axis = M * average/norm(average)
                    push!(axes_sym_plots, lines!(ax, [Point3f(-padding*rot_axis[:]), Point3f(padding*rot_axis[:])], color = (color_pool[1], dual_axes_opacity), linewidth = 4))
                end
            end

            points_plot = scatter!(ax, points, color = syms_color[1], markersize = 20)
            scatter!(ax, Point3f(0,0,0), color = movie_color_pool[1], markersize = 20, label = "Origin")
            lines!(ax, s_data[1,:], s_data[2,:], s_data[3,:], color = syms_color[1][1])

            # axislegend("Legend")
            
            been_primal = zeros(length(syms))

            function change_function(t)
                new_points = [Point3f(0,0,0) for _ in 1:num_gen]
                for i ∈ 1:num_gen
                    p = Point3f(gens[:,:,i] * s_data[:,t])

                    new_points[i] = p
                end

                points[] = new_points

                # Rank the symmetries
                ranking = rank_symmetries(new_points, syms)

                # Updating which symmetry to show
                for i ∈ eachindex(ranking)
                    sym_index = ranking[i]
                    sym = syms[sym_index]
                    shape_sides, num_shapes = size(sym)
                    # Primal
                    if i == 1
                        been_primal[sym_index] = 1
                        color_vec = Vector{}(undef, num_gen)
                        for j ∈ 1:num_gen
                            color_vec[sym[j]] = syms_color[sym_index][j]
                        end
                        points_plot.color = color_vec

                        for j ∈ 1:num_shapes
                            outline_plots[sym_index][j].color = (syms_color[sym_index][1 + (j-1)*shape_sides], primal_outline_opacity)
                        end

                        for j ∈ 1:num_shapes
                            mesh_plots[sym_index][j].color = (syms_color[sym_index][1 + (j-1)*shape_sides], primal_mesh_opacity)
                        end

                        for p ∈ axes_plots[sym_index]
                            p.color = (:black, primal_axes_opacity)
                        end
                    # Dual    
                    elseif i <= 1 + num_duals

                        for j ∈ 1:num_shapes
                            outline_plots[sym_index][j].color = (:black, dual_outline_opacity)
                        end

                        for j ∈ 1:num_shapes
                            mesh_plots[sym_index][j].color = (:black, dual_mesh_opacity)
                        end

                        for p ∈ axes_plots[sym_index]
                            p.color = (:black, dual_axes_opacity)
                        end
                    else

                        for j ∈ 1:num_shapes
                            outline_plots[sym_index][j].color = (:black, 0.0)
                        end

                        for j ∈ 1:num_shapes
                            mesh_plots[sym_index][j].color = (:black, 0.0)
                        end

                        for p ∈ axes_plots[sym_index]
                            p.color = (:black, 0.0)
                        end
                    end
                end

                # Determine current extent
                xs = getindex.(new_points, 1)
                ys = getindex.(new_points, 2)
                zs = getindex.(new_points, 3)

                if eccentric
                    M = maximum(abs, vcat(xs, ys, zs))

                    ax.limits = (
                        (-padding*M, padding*M),
                        (-padding*M, padding*M),
                        (-padding*M, padding*M)
                    )
                end

                # Save movie thumbnail
                if __save__


                end
            end

            if __save__
                record(change_function, fig, joinpath(__save_location__,  name * ".mp4"), time_step; framerate = framerate, px_per_unit = 2.0)
            else
                record(change_function, fig, joinpath(__save_location__,  "temp" * ".mp4"), time_step; framerate = framerate, px_per_unit = 2.0)
            end           
        end

        if picture

            fig = Figure(size = (800,800), padding = 0)
            ax = Axis3(fig[1,1], 
                perspectiveness = 0.5,
                azimuth = 1.3 * π,
                elevation = π/8,
                aspect = :data,
                xlabel = L"$x$",
                xlabelsize = 20,
                ylabel = L"$y$",
                ylabelsize = 20,
                zlabel = L"$z$",
                zlabelsize = 20,
                limits = ((-padding*M,padding*M), (-padding*M,padding*M), (-padding*M,padding*M))
            )

            # Generator
            gen = rot * (1+e) * [-1;0;0]

            points = Vector{Point3f}(undef, num_gen)

            for i ∈ 1:num_gen
                points[i] = gens[:,:,i] * gen
            end

            ranking = rank_symmetries(points, syms)
            r = 0.30*M/sqrt(num_gen)

            syms_color = []

            for i ∈ eachindex(ranking)
                sym_index = ranking[i]
                sym = syms[sym_index]
                shape_sides, num_shapes = size(sym)

                sym_colors = repeat(picture_color_pool[3:num_shapes+2], inner=shape_sides)
                push!(syms_color, sym_colors)

                # Primal
                if i == 1
                    plot_cloud(ax, points, sym, sym_colors, 1, r)
                    plot_polygons(ax, points, sym, sym_colors, primal_outline_opacity, primal_mesh_opacity)
                    plot_axes(ax, points, sym, primal_axes_opacity, padding)
                elseif i <= num_duals + 1
                    plot_polygons(ax, points, sym, sym_colors, dual_outline_opacity, dual_mesh_opacity)
                    plot_axes(ax, points, sym, dual_axes_opacity, padding)
                end
            end

            # Plotting orbit and generating_orbit
            meshscatter!(ax, [Point3f(0,0,0)], color = picture_color_pool[1], markersize = r)
            lines!(ax, generating_orbit[1,:], generating_orbit[2,:], generating_orbit[3,:], color = syms_color[1][1])

            buf = colorbuffer(fig)
            cropped_buf = buf[1 + Int(0.15 * size(buf,1)):end,:]

            if __display__
                display(fig)
            end

            if __save__
                save(joinpath(__save_location__, name * ".png"), cropped_buf)
            end

        end
    end
end

function triangulation(shape_sides)
    if shape_sides == 3
        return [
            GeometryBasics.TriangleFace(1, 2, 3)
        ]
    elseif shape_sides == 4
        return [
            GeometryBasics.TriangleFace(1, 2, 3),
            GeometryBasics.TriangleFace(1, 3, 4)
        ]
    elseif shape_sides == 5
        return [
            GeometryBasics.TriangleFace(1, 2, 3),
            GeometryBasics.TriangleFace(1, 3, 4),
            GeometryBasics.TriangleFace(1, 4, 5)
        ]
    end
end

function shade(c, centre; ambient = 0.38)
    LIGHT = normalize([0.6, 0.5, 0.8])
    t = ambient + (1-ambient)*max(0.0, dot(normalize(centre), LIGHT))
    # Convert named colors like :red back to their rgb values
    if typeof(c) == Symbol
        c = parse(RGB, String(c))
    end
    RGBf(t*red(c), t*green(c), t*blue(c))
end

function rank_symmetries(points, syms)
    sym_distances = zeros(length(syms))
    for i ∈ eachindex(syms)
        sym = syms[i]
        sym_distances[i] = norm(points[sym[1]] - points[sym[2]])
    end
    return sortperm(sym_distances)
end

function plot_cloud(ax, points, sym, colors, opacity, r)
    # Cloud of points
    color_vec = Vector{}(undef, length(points))
    for j ∈ eachindex(points)
        color_vec[sym[j]] = colors[j]
    end
    point_cloud = meshscatter!(ax, points, color = color_vec, markersize = r)
    point_cloud.opacity = opacity
    return point_cloud
end

function plot_polygons(ax, points, sym, colors, outline_opacity, mesh_opacity)
    # Polygon outlines and mesh
    shape_sides, num_shapes = size(sym)

    outline_plots = []
    mesh_plots = []
    for j = 1:num_shapes
        starting_index = 1 + (j-1)*shape_sides
        polygon_points = Vector{Point3f}(undef, shape_sides+1)
        for k = 1:shape_sides
            current_index = starting_index + k - 1
            polygon_points[k] = points[sym[current_index]]
        end
        polygon_points[shape_sides+1] = points[sym[starting_index]]
        push!(outline_plots, lines!(ax, polygon_points, color = (colors[starting_index], outline_opacity)))

        vertices = polygon_points[1:shape_sides]

        faces = triangulation(shape_sides)

        centre = vec(sum(vertices)) / shape_sides
        push!(mesh_plots, mesh!(ax, vertices, faces,
        color = (shade(colors[starting_index], centre), mesh_opacity)))
    end

    return outline_plots, mesh_plots
end

function plot_axes(ax, points, sym, opacity, padding)
    shape_sides, num_shapes = size(sym)
    for i = 1:num_shapes
        starting_index = 1 + (i-1)*shape_sides
        vertices = Vector{Point3f}(undef, shape_sides)
        for j = 1:shape_sides
            current_index = starting_index + j - 1
            vertices[j] = points[sym[current_index]]
        end
        average = 1 / shape_sides * sum(vertices)
        rot_ax = maximum(norm.(points)) * average / norm(average)

        lines!(ax, [-padding*rot_ax, padding*rot_ax], color = (:black, opacity), linewidth = 4)
    end
end


main()