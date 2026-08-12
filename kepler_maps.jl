function f_kepler(θ, c, e)
    return π * c^(-3) * (1 + e * cos(θ))^(-2)
end

function df_kepler(θ, c, e)
    return [
        2 * π * c^(-3) * (1 + e* cos(θ))^(-3) * sin(θ)
        -3 * π * c^(-4) * (1 + e * cos(θ))^(-2)
        ]
end

function K!(K, X, e)
    # Extract components of X
    θ = block(X,1)
    c = block(X,2)[1]

    Ncheb = order(θ)

    # Reset the container
    K.coefficients .= 0

    # Compute the Chebyshev coefficients of the nonlinearity
    G = zeros(space(θ))
    G_of_u!(G, θ, θ -> kepler_ode_map(θ, c, e), Ncheb)

    # Building diagonal operator D and shift operator T
    T = zeros(space(θ), space(θ))
    D = zeros(Chebyshev(Ncheb+1), space(θ))

    for k = 1:Ncheb
        D[k,k] = 2 * K
        T[k, k-1] = 1
        T[k, k+1] = -1
    end
    
    block(K,1).coefficients[:] = (D * θ - T * G).coefficients[:] 

    block(K,1)[0] = θ(-1)

    block(K,2)[1] = θ(1)


end

function DK!(DK, X, e, Ncheb)

end