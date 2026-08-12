# General code that will produce finite dimensional proofs

function first_order_radii_polymnomial(X, f::Function, Df::Function; r_star = 1e-5)
    iX = interval.(X)
    ir = interval(r_star)

    println("f(x) = ")
    display(mid.(f(iX)))

    A = interval.(inv(mid.(Df(iX))))

    println("A = ")
    display(A)

    F = f(iX)
    DF = Df(iX + ir*interval(-1,1)*one.(iX))

    Y = norm(A*F, Inf)
    Z = norm(interval(I) - A * DF, Inf)
    r = Y / (interval(1) - Z)

    println("\nY = " * string(Y))
    println("Z = " * string(Z))
    println("radius = " * string(sup(r)))

    if sup(Z) >= 1 || sup(r) >= inf(r_star)
        return (sup(r), inf(r_star), false)
    else
        return (sup(r), inf(r_star), true)
    end
end