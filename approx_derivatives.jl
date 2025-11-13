function DF_approx!(DF_approx::LinearOperator, u::Sequence, N_fft::Int64)
    S = space(component(u,1))
    N = order(S)
    h = 1e-6

    new_u = zeros(ComplexF64, S^3)

    F = zeros(ComplexF64, S^3)
    Fₕ1 = zeros(ComplexF64, S^3)
    Fₕ2 = zeros(ComplexF64, S^3)

    F!(F, u, N_fft)

    for i in 1:3
        for n = -N:N
            new_u[:] = u[:]
            component(new_u,i)[n] = component(u,i)[n] - h
            F!(Fₕ1, new_u, N_fft)

            component(new_u,i)[n] = component(u,i)[n] + h
            F!(Fₕ2, new_u, N_fft)

            approx = 1/(2*h) * (Fₕ2-Fₕ1);

            component(DF_approx,:,i)[:,n] = approx[:]
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

function DA_1_approx(ψ, e)
    h = 1e-6;

    approx = zeros(1,4)

    approx[i] = 1/2h * ()

end