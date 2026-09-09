include("action.jl")
include("proofs.jl")
include("algorithms.jl")

using LinearAlgebra
using GLMakie
using Printf


# ============================================================
# CANONICAL CYCLE
# ============================================================

function canonical_cycle(path)

    n = length(path)

    # Forward orientation
    min_position = argmin(path)

    forward = [
        path[
            mod1(
                min_position + k,
                n
            )
        ]
        for k in 0:n-1
    ]

    # Reverse orientation
    reverse_path = reverse(path)

    min_position_reverse =
        argmin(reverse_path)

    backward = [
        reverse_path[
            mod1(
                min_position_reverse + k,
                n
            )
        ]
        for k in 0:n-1
    ]

    if Tuple(forward) < Tuple(backward)
        return Tuple(forward)
    else
        return Tuple(backward)
    end

end


# ============================================================
# BUILD ADJACENCY GRAPH
#
# Two vertices are connected if their distance equals
# edge_length within tolerance.
# ============================================================

function build_adjacency(
    points,
    edge_length,
    tol
)

    num_vertices =
        size(points, 2)

    adjacency = [
        Int[]
        for _ in 1:num_vertices
    ]

    allowed_error =
        tol *
        max(
            1.0,
            abs(edge_length)
        )

    for i in 1:num_vertices

        for j in i+1:num_vertices

            d =
                norm(
                    points[:, i] -
                    points[:, j]
                )

            if abs(
                d - edge_length
            ) <= allowed_error

                push!(
                    adjacency[i],
                    j
                )

                push!(
                    adjacency[j],
                    i
                )

            end

        end

    end

    for i in 1:num_vertices
        sort!(adjacency[i])
    end

    return adjacency

end


# ============================================================
# FIND ALL SIMPLE CYCLES OF A GIVEN ORDER
#
# Searches ALL cycles of exactly cycle_size.
#
# Every edge belongs to the same distance class because the
# adjacency graph was constructed using one edge length.
# ============================================================

function find_cycles(
    adjacency,
    cycle_size
)

    num_vertices =
        length(adjacency)

    cycles =
        Set{
            Tuple{Vararg{Int}}
        }()


    function extend_cycle(path)

        start =
            path[1]

        current =
            path[end]


        # ----------------------------------------------------
        # Desired cycle order reached.
        # ----------------------------------------------------

        if length(path) ==
           cycle_size

            if start in adjacency[current]

                cycle =
                    canonical_cycle(path)

                push!(
                    cycles,
                    cycle
                )

            end

            return

        end


        # ----------------------------------------------------
        # Continue through every possible neighbour.
        # ----------------------------------------------------

        for next_vertex in
            adjacency[current]


            # Cannot reuse a vertex.
            if next_vertex in path
                continue
            end


            # Symmetry breaking.
            if next_vertex < start
                continue
            end


            push!(
                path,
                next_vertex
            )

            extend_cycle(path)

            pop!(
                path
            )

        end

    end


    # Start from every vertex.
    for start in 1:num_vertices

        extend_cycle(
            [start]
        )

    end


    return sort(
        collect(cycles)
    )

end


# ============================================================
# VERIFY CYCLE
# ============================================================

function verify_cycle(
    cycle,
    points,
    expected_length,
    tol
)

    n =
        length(cycle)

    allowed_error =
        tol *
        max(
            1.0,
            abs(expected_length)
        )

    for j in 1:n

        a =
            cycle[j]

        b =
            cycle[
                mod1(
                    j + 1,
                    n
                )
            ]

        d =
            norm(
                points[:, a] -
                points[:, b]
            )

        if abs(
            d - expected_length
        ) > allowed_error

            return false

        end

    end

    return true

end


# ============================================================
# WRITE SUMMARY CSV
#
# One row per distance class.
# Columns are cycle orders 3 through 20.
# ============================================================

function write_summary_csv(
    filename,
    distance_classes,
    adjacency_by_distance,
    all_cycles_by_distance,
    min_cycle_size,
    max_cycle_size
)

    open(
        filename,
        "w"
    ) do io

        # Header

        print(
            io,
            "distance_class,edge_length,number_of_edges"
        )

        for n in
            min_cycle_size:max_cycle_size

            print(
                io,
                ",cycles_",
                n
            )

        end

        println(io)


        # Data rows

        for d in
            1:length(distance_classes)


            edge_length =
                distance_classes[d]

            adjacency =
                adjacency_by_distance[d]

            cycles =
                all_cycles_by_distance[d]


            number_of_edges =
                sum(
                    length.(adjacency)
                ) ÷ 2


            print(
                io,
                d,
                ",",
                edge_length,
                ",",
                number_of_edges
            )


            for n in
                min_cycle_size:max_cycle_size

                print(
                    io,
                    ",",
                    length(
                        cycles[n]
                    )
                )

            end


            println(io)

        end

    end

end


# ============================================================
# WRITE NONZERO RESULTS CSV
#
# One row for every distance/order combination that actually
# contains at least one cycle.
# ============================================================

function write_nonzero_csv(
    filename,
    distance_classes,
    all_cycles_by_distance,
    min_cycle_size,
    max_cycle_size
)

    open(
        filename,
        "w"
    ) do io

        println(
            io,
            "distance_class,edge_length,cycle_order,number_of_cycles"
        )


        for d in
            1:length(distance_classes)


            edge_length =
                distance_classes[d]

            cycles =
                all_cycles_by_distance[d]


            for n in
                min_cycle_size:max_cycle_size


                number =
                    length(
                        cycles[n]
                    )


                if number > 0

                    println(
                        io,
                        d,
                        ",",
                        edge_length,
                        ",",
                        n,
                        ",",
                        number
                    )

                end

            end

        end

    end

end


# ============================================================
# WRITE INDIVIDUAL CYCLE CSV
#
# Every actual cycle gets one row.
# ============================================================

function write_cycles_csv(
    filename,
    distance_classes,
    all_cycles_by_distance,
    min_cycle_size,
    max_cycle_size
)

    open(
        filename,
        "w"
    ) do io

        println(
            io,
            "distance_class,edge_length,cycle_order,cycle_number,vertices"
        )


        for d in
            1:length(distance_classes)


            edge_length =
                distance_classes[d]

            cycles =
                all_cycles_by_distance[d]


            for n in
                min_cycle_size:max_cycle_size


                cycle_list =
                    cycles[n]


                for k in
                    1:length(cycle_list)


                    cycle =
                        cycle_list[k]


                    vertices =
                        join(
                            cycle,
                            ";"
                        )


                    println(
                        io,
                        d,
                        ",",
                        edge_length,
                        ",",
                        n,
                        ",",
                        k,
                        ",\"",
                        vertices,
                        "\""
                    )

                end

            end

        end

    end

end


# ============================================================
# PRINT COMPACT SUMMARY TABLE
# ============================================================

function print_summary_table(
    distance_classes,
    adjacency_by_distance,
    all_cycles_by_distance,
    min_cycle_size,
    max_cycle_size
)

    println()
    println("============================================================")
    println("DISTANCE × CYCLE-ORDER SUMMARY")
    println("============================================================")

    println()
    println(
        "Each row represents one distinct edge length."
    )

    println(
        "Each entry is the number of simple cycles of that order."
    )

    println()


    # --------------------------------------------------------
    # Header
    # --------------------------------------------------------

    @printf(
        "%8s %16s %8s",
        "Class",
        "Distance",
        "Edges"
    )

    for n in
        min_cycle_size:max_cycle_size

        @printf(
            " %8d",
            n
        )

    end

    println()


    # Separator

    width =
        8 +
        1 +
        16 +
        1 +
        8 +
        1 +
        8 *
        (
            max_cycle_size -
            min_cycle_size +
            1
        )

    println(
        "-"^width
    )


    # --------------------------------------------------------
    # Rows
    # --------------------------------------------------------

    for d in
        1:length(distance_classes)


        edge_length =
            distance_classes[d]

        adjacency =
            adjacency_by_distance[d]

        cycles =
            all_cycles_by_distance[d]


        number_of_edges =
            sum(
                length.(adjacency)
            ) ÷ 2


        @printf(
            "%8d %16.10f %8d",
            d,
            edge_length,
            number_of_edges
        )


        for n in
            min_cycle_size:max_cycle_size


            @printf(
                " %8d",
                length(
                    cycles[n]
                )
            )

        end


        println()

    end

end


# ============================================================
# PRINT ONLY DISTANCES THAT CONTAIN CYCLES
# ============================================================

function print_nonzero_table(
    distance_classes,
    all_cycles_by_distance,
    min_cycle_size,
    max_cycle_size
)

    println()
    println()
    println("============================================================")
    println("DISTANCES CONTAINING CYCLES")
    println("============================================================")

    println()
    println(
        "Only nonzero cycle counts are shown."
    )

    println()


    for d in
        1:length(distance_classes)


        edge_length =
            distance_classes[d]

        cycles =
            all_cycles_by_distance[d]


        found_something =
            false


        for n in
            min_cycle_size:max_cycle_size


            number =
                length(
                    cycles[n]
                )


            if number > 0


                if !found_something

                    println()
                    println(
                        "Distance class ",
                        d,
                        "   distance = ",
                        @sprintf(
                            "%.10f",
                            edge_length
                        )
                    )

                    found_something =
                        true

                end


                println(
                    "    order ",
                    n,
                    " : ",
                    number,
                    " cycle(s)"
                )

            end

        end

    end

end


# ============================================================
# PRINT TOTAL NUMBER OF CYCLES BY ORDER
# ============================================================

function print_totals_by_order(
    distance_classes,
    all_cycles_by_distance,
    min_cycle_size,
    max_cycle_size
)

    println()
    println()
    println("============================================================")
    println("TOTAL CYCLES BY ORDER")
    println("============================================================")

    println()


    for n in
        min_cycle_size:max_cycle_size


        total =
            sum(
                length(
                    all_cycles_by_distance[d][n]
                )
                for d in
                1:length(distance_classes)
            )


        @printf(
            "Order %2d : %10d cycle(s)\n",
            n,
            total
        )

    end

end


# ============================================================
# PRINT WHICH DISTANCES PRODUCE EACH ORDER
# ============================================================

function print_orders_by_distance(
    distance_classes,
    all_cycles_by_distance,
    min_cycle_size,
    max_cycle_size
)

    println()
    println()
    println("============================================================")
    println("DISTANCES PRODUCING EACH CYCLE ORDER")
    println("============================================================")

    println()


    for n in
        min_cycle_size:max_cycle_size


        entries =
            Tuple{Int,Float64,Int}[]


        for d in
            1:length(distance_classes)


            number =
                length(
                    all_cycles_by_distance[d][n]
                )


            if number > 0

                push!(
                    entries,
                    (
                        d,
                        distance_classes[d],
                        number
                    )
                )

            end

        end


        if !isempty(entries)


            println(
                "ORDER ",
                n
            )

            println(
                "------------------------------------------------------------"
            )


            for (
                d,
                distance,
                number
            ) in entries


                println(
                    "    class ",
                    d,
                    "   distance = ",
                    @sprintf(
                        "%.10f",
                        distance
                    ),
                    "   cycles = ",
                    number
                )

            end


            println()

        end

    end

end


# ============================================================
# PRINT ACTUAL CYCLES FOR A SPECIFIC DISTANCE AND ORDER
#
# Example:
#
# print_cycles(
#     results,
#     7,
#     6
# )
#
# prints all 6-cycles at distance class 7.
# ============================================================

function print_cycles(
    results,
    distance_class,
    cycle_order
)

    distance_classes =
        results.distance_classes

    all_cycles_by_distance =
        results.all_cycles_by_distance


    if distance_class < 1 ||
       distance_class > length(distance_classes)

        error(
            "Invalid distance class."
        )

    end


    if cycle_order < 3 ||
       cycle_order > 20

        error(
            "Cycle order must be between 3 and 20."
        )

    end


    edge_length =
        distance_classes[
            distance_class
        ]


    cycles =
        all_cycles_by_distance[
            distance_class
        ][cycle_order]


    println()
    println("============================================================")
    println("CYCLES")
    println("============================================================")

    println(
        "Distance class = ",
        distance_class
    )

    println(
        "Edge length = ",
        edge_length
    )

    println(
        "Cycle order = ",
        cycle_order
    )

    println(
        "Number of cycles = ",
        length(cycles)
    )

    println()


    for k in
        1:length(cycles)

        println(
            "Cycle ",
            k,
            " = ",
            collect(
                cycles[k]
            )
        )

    end

end


# ============================================================
# MAIN SEARCH
# ============================================================

function run_cycle_search()


    # ========================================================
    # GENERATORS
    # ========================================================

    gens, num_gen =
        generators_60()


    # ========================================================
    # ORBIT PARAMETERS
    # ========================================================

    X = [
        5.793233619038003; 6.190758148847939; 6.25975674713244; 0.6795057169690465
    ]

    ψ =
        X[1:3]

    e =
        X[4]

    rot =
        rotation_matrix(ψ)

    θ = 0


    # ========================================================
    # INITIAL VECTOR
    # ========================================================

    v =
        rot *
        (1 - e^2) /
        (1 + e*cos(θ)) *
        [cos(θ), sin(θ), 0]


    println()
    println("============================================================")
    println("ORBIT")
    println("============================================================")

    println(
        "Number of points = ",
        num_gen
    )

    println(
        "Cycle orders searched = 3 through ",
        min(
            20,
            num_gen
        )
    )

    println()


    display(v)


    # ========================================================
    # GENERATE POINTS
    # ========================================================

    points =
        zeros(
            3,
            num_gen
        )


    for i in
        1:num_gen

        points[:, i] =
            gens[:, :, i] *
            v

    end


    display(
        GLMakie.scatter(points)
    )


    # ========================================================
    # SETTINGS
    # ========================================================

    tol =
        1e-7

    min_cycle_size =
        3

    max_cycle_size =
        min(
            20,
            num_gen
        )


    # ========================================================
    # CALCULATE ALL PAIRWISE DISTANCES
    # ========================================================

    all_distances =
        Float64[]


    for i in
        1:num_gen

        for j in
            i+1:num_gen


            d =
                norm(
                    points[:, i] -
                    points[:, j]
                )


            push!(
                all_distances,
                d
            )

        end

    end


    sort!(
        all_distances
    )


    # ========================================================
    # GROUP EQUAL DISTANCES
    # ========================================================

    distance_classes =
        Float64[]


    for d in
        all_distances


        if isempty(
            distance_classes
        )


            push!(
                distance_classes,
                d
            )


        else


            reference =
                distance_classes[end]


            allowed_error =
                tol *
                max(
                    1.0,
                    abs(reference)
                )


            if abs(
                d - reference
            ) > allowed_error


                push!(
                    distance_classes,
                    d
                )

            end

        end

    end


    # ========================================================
    # PRINT DISTANCE CLASSES
    # ========================================================

    println()
    println("============================================================")
    println("DISTANCE CLASSES")
    println("============================================================")

    println(
        "Found ",
        length(distance_classes),
        " distinct distances."
    )

    println()


    for d in
        1:length(distance_classes)


        println(
            lpad(
                string(d),
                5
            ),
            " : ",
            @sprintf(
                "%.12f",
                distance_classes[d]
            )
        )

    end


    # ========================================================
    # STORAGE
    # ========================================================

    adjacency_by_distance =
        Dict()

    all_cycles_by_distance =
        Dict()


    # ========================================================
    # SEARCH EVERY DISTANCE
    #
    # AND EVERY ORDER 3 THROUGH 20.
    # ========================================================

    println()
    println("============================================================")
    println("SEARCHING")
    println("============================================================")


    for d in
        1:length(distance_classes)


        edge_length =
            distance_classes[d]


        println()
        println(
            "Distance class ",
            d,
            " / ",
            length(distance_classes),
            "   distance = ",
            @sprintf(
                "%.10f",
                edge_length
            )
        )


        # ----------------------------------------------------
        # Build equal-length graph.
        # ----------------------------------------------------

        adjacency =
            build_adjacency(
                points,
                edge_length,
                tol
            )


        adjacency_by_distance[d] =
            adjacency


        number_of_edges =
            sum(
                length.(adjacency)
            ) ÷ 2


        println(
            "    edges = ",
            number_of_edges
        )


        # ----------------------------------------------------
        # Search EVERY order.
        # ----------------------------------------------------

        cycles_for_distance =
            Dict()


        for n in
            min_cycle_size:max_cycle_size


            cycles =
                find_cycles(
                    adjacency,
                    n
                )


            cycles_for_distance[n] =
                cycles


            if !isempty(cycles)

                println(
                    "    order ",
                    n,
                    " = ",
                    length(cycles),
                    " cycle(s)"
                )

            end

        end


        all_cycles_by_distance[d] =
            cycles_for_distance

    end


    # ========================================================
    # SUMMARY TABLES
    # ========================================================

    print_summary_table(
        distance_classes,
        adjacency_by_distance,
        all_cycles_by_distance,
        min_cycle_size,
        max_cycle_size
    )


    print_nonzero_table(
        distance_classes,
        all_cycles_by_distance,
        min_cycle_size,
        max_cycle_size
    )


    print_totals_by_order(
        distance_classes,
        all_cycles_by_distance,
        min_cycle_size,
        max_cycle_size
    )


    print_orders_by_distance(
        distance_classes,
        all_cycles_by_distance,
        min_cycle_size,
        max_cycle_size
    )


    # ========================================================
    # VERIFY ALL CYCLES
    # ========================================================

    println()
    println()
    println("============================================================")
    println("VERIFYING CYCLES")
    println("============================================================")


    total_cycles =
        0

    verification_failed =
        false


    for d in
        1:length(distance_classes)


        edge_length =
            distance_classes[d]


        cycles_for_distance =
            all_cycles_by_distance[d]


        for n in
            min_cycle_size:max_cycle_size


            cycles =
                cycles_for_distance[n]


            for cycle in
                cycles


                total_cycles +=
                    1


                valid =
                    verify_cycle(
                        cycle,
                        points,
                        edge_length,
                        tol
                    )


                if !valid


                    verification_failed =
                        true


                    println(
                        "FAILED: distance class ",
                        d,
                        ", order ",
                        n,
                        ", cycle ",
                        collect(cycle)
                    )

                end

            end

        end

    end


    if verification_failed

        println()
        println(
            "WARNING: verification failed for some cycles."
        )

    else

        println()
        println(
            "✓ All ",
            total_cycles,
            " cycles passed verification."
        )

    end


    # ========================================================
    # WRITE CSV FILES
    # ========================================================

    summary_filename =
        "I_cycle_summary.csv"

    nonzero_filename =
        "I_cycle_nonzero.csv"

    cycles_filename =
        "I_cycles_all.csv"


    write_summary_csv(
        summary_filename,
        distance_classes,
        adjacency_by_distance,
        all_cycles_by_distance,
        min_cycle_size,
        max_cycle_size
    )


    write_nonzero_csv(
        nonzero_filename,
        distance_classes,
        all_cycles_by_distance,
        min_cycle_size,
        max_cycle_size
    )


    write_cycles_csv(
        cycles_filename,
        distance_classes,
        all_cycles_by_distance,
        min_cycle_size,
        max_cycle_size
    )


    # ========================================================
    # FINAL OUTPUT
    # ========================================================

    println()
    println()
    println("============================================================")
    println("FILES WRITTEN")
    println("============================================================")

    println()
    println(
        "1. ",
        summary_filename
    )

    println(
        "   Complete distance × order table."
    )

    println()
    println(
        "2. ",
        nonzero_filename
    )

    println(
        "   Only distance/order combinations containing cycles."
    )

    println()
    println(
        "3. ",
        cycles_filename
    )

    println(
        "   Every individual cycle and its vertex indices."
    )


    # ========================================================
    # FINAL SUMMARY
    # ========================================================

    println()
    println("============================================================")
    println("FINAL SUMMARY")
    println("============================================================")

    println(
        "Points             : ",
        num_gen
    )

    println(
        "Distance classes   : ",
        length(distance_classes)
    )

    println(
        "Cycle orders       : ",
        min_cycle_size,
        " - ",
        max_cycle_size
    )

    println(
        "Total cycles       : ",
        total_cycles
    )

    println()
    println(
        "Search complete."
    )


    # ========================================================
    # RETURN RESULTS
    # ========================================================

    return (
        points = points,

        distance_classes =
            distance_classes,

        adjacency_by_distance =
            adjacency_by_distance,

        all_cycles_by_distance =
            all_cycles_by_distance
    )

end


# ============================================================
# RUN
# ============================================================

results =
    run_cycle_search()


# ============================================================
# CONVENIENCE VARIABLES
# ============================================================

points =
    results.points

distance_classes =
    results.distance_classes

adjacency_by_distance =
    results.adjacency_by_distance

all_cycles_by_distance =
    results.all_cycles_by_distance


# ============================================================
# EXAMPLES OF HOW TO INSPECT SPECIFIC RESULTS
# ============================================================

# To see all 3-cycles at distance class 5:
#
# print_cycles(
#     results,
#     5,
#     3
# )


# To see all 7-cycles at distance class 3:
#
# print_cycles(
#     results,
#     3,
#     7
# )


# To access the cycles directly:
#
# all_cycles_by_distance[5][3]
#
# gives all 3-cycles at distance class 5.
#
# all_cycles_by_distance[3][7]
#
# gives all 7-cycles at distance class 3.
