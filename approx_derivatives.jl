function DF_approx!(DF_approx::LinearOperator, u::Sequence, ε::Float64, N_fft::Int64)
    S = space(component(u,1))
    N = order(S)
    h = 1e-6

    new_u = zeros(ComplexF64, S^3)

    F = zeros(ComplexF64, S^3)
    Fₕ1 = zeros(ComplexF64, S^3)
    Fₕ2 = zeros(ComplexF64, S^3)

    F!(F, u, ε, N_fft)

    for i in 1:3
        for n = -N:N
            new_u[:] = u[:]
            component(new_u,i)[n] = component(u,i)[n] - h
            F!(Fₕ1, new_u, ε, N_fft)

            component(new_u,i)[n] = component(u,i)[n] + h
            F!(Fₕ2, new_u, ε, N_fft)

            approx = 1/(2*h) * (Fₕ2-Fₕ1);

            component(DF_approx,:,i)[:,n] = approx[:]
        end
    end
end

function Dg_approx(u)
    h = 1e-6;

    approx = zeros(3,3)

    approx[:,1] = 1/(2*h) * (g(u + h*[1;0;0])-g(u - h*[1;0;0]))
    approx[:,2] = 1/(2*h) * (g(u + h*[0;1;0])-g(u - h*[0;1;0]))
    approx[:,3] = 1/(2*h) * (g(u + h*[0;0;1])-g(u - h*[0;0;1]))

    return approx
end

function Dh_approx(u)
    k = 1e-6;

    approx = zeros(3,3)

    approx[:,1] = 1/(2*k) * (h(u + k*[1;0;0])-h(u - k*[1;0;0]))
    approx[:,2] = 1/(2*k) * (h(u + k*[0;1;0])-h(u - k*[0;1;0]))
    approx[:,3] = 1/(2*k) * (h(u + k*[0;0;1])-h(u - k*[0;0;1]))

    return approx
end

function DA_1_approx(X, N_fft)
    h = 1e-6;

    approx = zeros(4,1)

    approx[1] = 1/2h * (A_1(X + h*[1;0;0;0], N_fft) - A_1(X - h*[1;0;0;0], N_fft))
    approx[2] = 1/2h * (A_1(X + h*[0;1;0;0], N_fft) - A_1(X - h*[0;1;0;0], N_fft))
    approx[3] = 1/2h * (A_1(X + h*[0;0;1;0], N_fft) - A_1(X - h*[0;0;1;0], N_fft))
    approx[4] = 1/2h * (A_1(X + h*[0;0;0;1], N_fft) - A_1(X - h*[0;0;0;1], N_fft))

    return approx
end

function HA_1_approx(X, N_fft)
    h = 1e-6;

    approx = zeros(4,4)

    approx[:,1] = 1/2h * (DA_1(X + h*[1;0;0;0], N_fft) - DA_1(X - h*[1;0;0;0], N_fft))
    approx[:,2] = 1/2h * (DA_1(X + h*[0;1;0;0], N_fft) - DA_1(X - h*[0;1;0;0], N_fft))
    approx[:,3] = 1/2h * (DA_1(X + h*[0;0;1;0], N_fft) - DA_1(X - h*[0;0;1;0], N_fft))
    approx[:,4] = 1/2h * (DA_1(X + h*[0;0;0;1], N_fft) - DA_1(X - h*[0;0;0;1], N_fft))

    return approx
end