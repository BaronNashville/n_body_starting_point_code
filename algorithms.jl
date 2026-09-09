function MyNewton!(F_DF!, x, F, DF; tol=1e-13, max_iter=20)
    F_DF!(F, DF, x)
    i = 0

    println("Beginning Newton")
    println("Iteration " * string(i) * ", ||F(u)|| = " * string(norm(F)) * ", ||DF\\F|| = " * string(norm(DF \ F)))
    while i <= max_iter - 1 && norm(F) > tol
        i = i + 1
        x = x - DF \ F
        F_DF!(F, DF, x)
        println("Iteration " * string(i) * ", ||F(u)|| = " * string(norm(F)) * ", ||DF\\F|| = " * string(norm(DF \ F)))
    end
    println("Newton ended after " * string(i) * " iterations. \n")
    return x
end

function MyNewton(F_DF, x; tol=1e-13, max_iter=20, verbose = true)
    # We expect function input to be as follows
    # F_DF = x -> (F(x), DF(x))
    if verbose
        println("Beginning Newton")
    end

    Fx, DFx = F_DF(x)

    if ! ( iszero(isnan.(Fx)) && iszero(isnan.(DFx)) && iszero(isinf.(Fx)) && iszero(isinf.(DFx)))
        println("Newton stopped due to NaN or Inf values")
        return x, false
    end
    i = 0
    if verbose
        println("Iteration " * string(i) * ", ||F(u)|| = " * string(norm(Fx)) * ", u = " * string(x))
    end
    while i <= max_iter - 1 && norm(Fx) > tol
        i = i + 1
        try
            x = (x-DFx\Fx)[:]
        catch e
            println("Newton stopped due to error. Error: " * string(e))
            return x, false
        end

        Fx, DFx = F_DF(x)

        if !( iszero(isnan.(Fx)) && iszero(isnan.(DFx)) && iszero(isinf.(Fx)) && iszero(isinf.(DFx)))
            println("Newton stopped due to NaN or Inf values")
            return x, false
        end
        if verbose
            println("Iteration " * string(i) * ", ||F(u)|| = " * string(norm(Fx)) * ", ||DF\\F|| = " * string(norm(DFx \ Fx)) * ", u = " * string(x))
        end
    end
    if verbose
        println("Newton ended after " * string(i) * " iterations. \n")
    end
    if norm(Fx) <= tol
        return x, true
    else
        return x, false
    end
end

function IterativeGradientDescent(F_DF, x; max_step_size::Float64=0.1, tol::Float64=1e-12, max_iter::Int64=1000, verbose = true)
    # We expect function input to be as follows
    # F_DF = x -> (F(x), DF(x))
    step_size = max_step_size
    count = 0
    Fx, DFx = F_DF(x)
    if verbose
        println("Beginning gradient descent")
    end
    while norm(DFx, Inf) > tol && count < max_iter
        step = (-step_size*DFx)[:]
        Fstep, DFstep = F_DF(x + step)
        # Checking qualify of test point
        if Fstep > Fx + 0.2 * step_size * norm(DFx, 2)
            step_size = step_size / 2
            continue
        end
        if verbose
            println("Iteration " * string(count) * ", u = " * string(x) * ", ||F(u)|| = " * string(norm(Fx)) * ", ||DF(u)|| = " * string(norm(DFx, Inf)))
        end
        x = x + step
        count = count + 1
        Fx, DFx = Fstep, DFstep
    end
    if verbose
        println("Iteration " * string(count) * ", u = " * string(x) * ", ||F(u)|| = " * string(norm(Fx)) * ", ||DF(u)|| = " * string(opnorm(DFx, Inf)) * "\nGradient descent ended after " * string(count) * " iterations. \n")
    end
    if norm(DFx, Inf) <= tol
        return x, true
    else
        return x, false
    end
    return x, success
end

function IntegralGradientDescent(DF, x; tol::Float64=1e-12, verbose = true)
    u₀ = x
    tspan = (0.0, 10.0)
    prob = ODEProblem((u, p, t) -> -DF(u)[:], u₀, tspan)

    solution = solve(prob, Tsit5())
    x = (solution.u)[end]

    if verbose
        println(solution.u)
    end

    if norm(DF(x), Inf) <= tol
        return x, true
    else
        return x, false
    end
end