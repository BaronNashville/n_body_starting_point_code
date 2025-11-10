function F!(F::Sequence, u::Sequence, N_fft::Int64)
    # Extract important space data
    S = space(component(u,1))
    G = zeros(ComplexF64, S^3)    
    N = order(S)

    # Initializing the output to be 0
    F .= 0

    # Linear part
    D = zeros(S^3,S^3)
    for i = 1:3
        for k in -N:N
            component(D,i,i)[k,k] = -k^2
        end
    end

    # Nonlinear part
    G_of_u_vec!(G, u, g, N_fft)

    F[:] = (D*u + G)[:] 
end

function DF!(DF::LinearOperator, u::Sequence, N_fft::Int64)
    # Extract important space data
    S = space(component(u,1))
    DG = zeros(ComplexF64, S^3,S^3)    
    N = order(S)

    # Initializing the output to be 0
    DF .= 0

    # Linear part
    D = zeros(S^3,S^3)
    for i = 1:3
        for k in -N:N
            component(D,i,i)[k,k] = -k^2
        end
    end

    # Nonlinear part
    G_of_u_mat!(DG, u, Dg, N_fft)

    DF.coefficients[:] = (D + DG).coefficients[:]  
end

function DG!(DG::LinearOperator, u::Sequence)
    # Extract important space data
    S = space(component(u,1))
    N = order(S)

    N_fft = NFFT(N)
    S_pad = Fourier(N_fft,frequency(S))

    # Project sequence into higher order fourier space by padding with zeros
    u_pad = project(u, S_pad^3)

    grid_points = zeros(ComplexF64, 3, 2*N_fft+1)
    for i = 1:3
        grid_points[i,:] = FFTW.fft(component(u_pad, i).coefficients[:])
    end

    eval_values = zeros(ComplexF64, 3, 3, 2*N_fft+1)
    for i = 1:2*N_fft+1
        eval_values[:,:,i] = Dg(grid_points[:,i])
    end

    for i = 1:3
        for j = 1:3
            component(DG,i,j).coefficients[:] = project(Multiplication(Sequence(S_pad, FFTW.ifft(eval_values[i,j,:]))), S, S).coefficients[:]
        end
    end

end

function g(u)
    return [exp(u[1]); exp(u[2]); exp(u[3])]
end

function Dg(u)
    return [exp(u[1]) 0 0; 0 exp(u[2]) 0; 0 0 exp(u[3])]
end

# function g(u)
#     return u/norm(u)^3
# end

# function Dg(u)
#     return 1/norm(u)^3 * I -3*u*transpose(u)/norm(u)^5
# end