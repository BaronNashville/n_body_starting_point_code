include("action.jl")
include("algorithms.jl")

N_fft  = 2^14
partitions = [3;3;3;5]
solutions = nothing

X = zeros(4)
for i ∈ 1:partitions[1]-1
    X[1] = 2 * π * (i / partitions[1])
    for j ∈ 1:partitions[2]-1
        X[2] = 2 * π * (j / partitions[2])
        for k ∈ 1:partitions[3]-1
            X[3] = 2 * π * (k / partitions[3])
            for ℓ ∈ 1:partitions[4]-1
                X[4] = ℓ / partitions[4]
                println("Current iteration: i = " * string(i) * ", j = " * string(j) * ", k = " * string(k) * ", ℓ = " * string(ℓ))

                # Newton approach
                sol, success = MyNewton(X -> (DA_1_float(X, N_fft, generators_24), HA_1_float(X, N_fft, generators_24)), X, tol = 1e-14, max_iter = 20, verbose = false)

                # Gradient descent method
                # sol, success = IterativeGradientDescent(X -> (A_1_float(X, N_fft), DA_1_float(X, N_fft)), X, tol = 1e-12, max_iter = 100, verbose = true)

                # Integrating ODE method
                # sol, success = IntegralGradientDescent(X -> DA_1_float(X, N_fft), X, tol = 1e-12, verbose = true)
                if success
                    global solutions
                    if abs(sol[end]) > 1e-4
                        if isnothing(solutions)
                            solutions = hcat(sol)
                        else
                            new = true
                            for i in size(solutions,2)
                                if norm(sol - solutions[:,i]) < 0.1
                                    new = false
                                    break
                                end
                            end

                            if new
                                solutions = hcat(solutions, sol)
                            end
                        end
                    end
                end
            end
        end
    end
end