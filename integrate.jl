# Integrates Qπ periodic function using the fft

############################################################
# Non rigorous implementation of integration using the fft #
############################################################
function fft_integrate_float(f::Function, Q, N_fft)
    # Create the grid
    period = Q * (π)
    # θ = zeros(N_fft)

    # for j = 0:N_fft-1
    #     θ[j+1] = period * big(j) / N_fft
    # end
    θ = (period / N_fft) * collect(0:N_fft-1)
    # Evaluate f at the grid
    f_θ = f.(θ)

    # Extract the 0th Fourier coefficient
    avg = (1) / N_fft * sum(f_θ)

    # Use it to compute the numerical integral
    value = avg * period

    return value
end

function fft_integrate(f::Function, Q::Rational{Int64}, N_fft::Int64; rigorous::Bool=false, ρ::Interval{Float64}=interval(0.01))
    # Create the grid
    # println("Creating the grid")
    if rigorous
        period = interval(Q) * interval(π)
        θ = interval.(zeros(N_fft))
    else
        period = Q * π
        θ = zeros(N_fft)
    end

    for j = 0:N_fft-1
        θ[j+1] = period * exact(j // N_fft)
    end
    # println("Grid\n" * string(θ))
    # Evaluate f at the grid
    # println("Sampling f at the grid")
    f_θ = f.(θ)
    # println("\nEvaluating f at the grid\n" * string(f_θ))
    # println("\nMax radius of f on the grid = " * string(maximum(radius.(f_θ))))

    # Extract the 0th Fourier coefficient
    # println("Computing average value")
    avg = exact(1 // N_fft) * sum(f_θ)
    # println("Radius of the average of f on the grid = " * string(radius(avg)))

    # Use it to compute the numerical integral
    # println("Use average to compute integral")
    value = avg * period
    # println("Numerical integral = " * string(value))
    # println("Radius of the integral = " * string(radius(value)))

    if rigorous
        # println("Validating the integral")
        # Verify the function has no poles in the required rectangle
        horizontal_grid = 100
        horizontal_length = interval(Q) * π / interval(horizontal_grid)
        vertical_grid = round(horizontal_grid * ρ / (Q*π))
        vertical_length = ρ / interval(vertical_grid)

        for i ∈ 0:horizontal_grid-1
            for j ∈ 0:vertical_grid-1
                box_up = horizontal_length * interval(i,i+1) + vertical_length * interval(1im) * interval(j, j+1) 
                box_down =  horizontal_length * interval(i,i+1) - vertical_length * interval(1im) * interval(j, j+1) 

                if sup(norm(f(box_up))) == Inf || sup(norm(f(box_down))) == Inf
                    println("Invalid choice of ρ, pole detected")
                end
            end
        end


        up_f = θ -> abs(f(θ + interval(1im) * ρ))
        down_f = θ -> abs(f(θ - interval(1im) * ρ))

        # println("Evaluating up shifted integral")
        C_up = lazy_integrate(up_f, interval(Q) * π, 1000, rigorous=true)
        # println("Evaluating down shifted integral")
        C_down = lazy_integrate(down_f, interval(Q) * π, 1000, rigorous=true)
        # println("Defining constant C")
        C = max(C_up, C_down)
        # println("C = " * string(C))

        # println("Defining constant nu")
        ν = exp(ρ)
        # println("Building the error")
        error = interval(sup(C * period * interval(2) / (ν^N_fft - interval(1))))
        # println("Aliassing error = " * string(error))
        value = value + error * interval(-1, 1) * period
    end

    return value
end

function lazy_integrate(f::Function, L::Real, N::Int64; rigorous::Bool=false)
    # Create the grid
    θ = zeros(typeof(L), N)
    for j = 0:N-1
        θ[j+1] = L * exact(j // N)
    end

    if rigorous
        # Turn grid into intervals covering the region of integration
        θ = θ .+ interval(0, L * exact(1 // N))
    end

    # Evaluate f at the grid
    f_θ = abs.(f.(θ))
    # println("\nEvaluating f at the grid\n" * string(maximum(f_θ)))

    return L * exact(1 // N) * sum(f_θ)
end