function DF_approx!(DF_approx::LinearOperator, u::Sequence, N_fft::Int64)
    S = space(component(u,1))
    N = order(S)
    h = 1e-6

    new_u = zeros(ComplexF64, S^3)

    F = zeros(ComplexF64, S^3)
    Fₕ = zeros(ComplexF64, S^3)

    F!(F, u, N_fft)

    for i in -N:N
        for j = 1:3
            new_u[:] = u[:]
            component(new_u,j)[i] = component(new_u,j)[i] + h

            F!(Fₕ, new_u, N_fft)

            approx = 1/h * (Fₕ-F);

            component(DF_approx,:,j)[:,i] = approx[:]
        end
    end
end

function Dg_approx(u)
    h = 1e-6;

    approx = zeros(3,3)

    approx[:,1] = 1/h * (g(u + h*[1;0;0])-g(u))
    approx[:,2] = 1/h * (g(u + h*[0;1;0])-g(u))
    approx[:,3] = 1/h * (g(u + h*[0;0;1])-g(u))

    return approx
end