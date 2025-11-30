function Newton!(F_DF!, x , F, DF; tol = 1e-12, max_iter = 20)
    F_DF!(F, DF, x)
    i = 0;

    println("Beginning Newton")
    println("Iteration " * string(i) * ", ||F(u)|| = " * string(norm(F)) * ", ||DF\\F|| = " * string(norm(DF \ F)))
    while i <= max_iter-1 && norm(F) > tol
        i = i+1
        x = x - DF \ F
        F_DF!(F, DF, x)
        println("Iteration " * string(i) * ", ||F(u)|| = " * string(norm(F)) * ", ||DF\\F|| = " * string(norm(DF \ F)))
    end
    println("Newton ended after " * string(i) * " iterations. \n")
    return x
end