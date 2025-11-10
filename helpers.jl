function G_of_u!(G_of_u::Sequence, u::Sequence, G::Function, N_fft::Int64)
    # INPUTS:
    # 1) u  = {u_n}_{n=-N}^{N} sequence of fourier coefficients
    # 2) g : R → R analytic function

    # Set ouput to 0
    G_of_u .= 0

    # Extract parameters from inputs
    S = space(u)
    f = frequency(S)

    S_pad = Fourier(div(N_fft,2), f)    

    u_pad = project(u, S_pad)

    # Evaluate the fourier series at grid points using the IFFT
    grid_points = Sequence(S_pad, [N_fft * FFTW.ifft(FFTW.ifftshift(u_pad.coefficients[begin:end-1])); 0])

    # Evaluate the function at the grid points
    eval_points = zeros(ComplexF64, S_pad)

    for n ∈ -div(N_fft,2):div(N_fft,2)-1
        eval_points[n] = G(grid_points[n])
    end

    # Use the FFT to obtain the coefficients of G(u)
    G_of_u.coefficients[:] = project(Sequence(S_pad, [1/N_fft * FFTW.fftshift(FFTW.fft(eval_points.coefficients[begin:end-1]));0]), S).coefficients[:]
end

function G_of_u_vec!(G_of_u::Sequence, u::Sequence, G::Function, N_fft::Int64)
    # INPUTS:
    # 1) u  = {u_n}_{n=-N}^{N} sequence of fourier coefficients
    # 2) g : Rⁿ → Rⁿ analytic function

    # Set ouput to 0
    G_of_u .= 0

    # Extract parameters from inputs
    S = space(component(u,1))
    f = frequency(S)
    dim = div(dimension(space(u)), dimension(S))

    S_pad = Fourier(div(N_fft,2), f)    

    u_pad = project(u, S_pad^dim)

    # Evaluate the fourier series at grid points using the IFFT
    grid_points = zeros(ComplexF64, dim, N_fft)
    for i ∈ 1:dim
        grid_points[i,:] = N_fft * FFTW.ifft(FFTW.ifftshift(component(u_pad,i).coefficients[begin:end-1]))
    end

    # Evaluate the function at the grid points
    eval_points = zeros(ComplexF64, dim, N_fft)

    for j ∈ 1:N_fft
        eval_points[:,j] = G(grid_points[:,j])
    end

    # Use the FFT to obtain the coefficients of G(u)
    fourier_coeffs = zeros(ComplexF64, dim, N_fft)
    for i ∈ 1:dim
        fourier_coeffs[i,:] = 1/N_fft * FFTW.fftshift(FFTW.fft(eval_points[i,:]))
    end

    for i ∈ 1:dim
        component(G_of_u,i).coefficients[:] = project(Sequence(S_pad, [fourier_coeffs[i,:];0]), S).coefficients[:]
    end  
end

function G_of_u_mat!(G_of_u::LinearOperator, u::Sequence, G::Function, N_fft::Int64)
    # INPUTS:
    # 1) u  = {u_n}_{n=-N}^{N} sequence of fourier coefficients
    # 2) G : Rⁿ → Rⁿ x Rⁿ  analytic function

    # Set ouput to 0
    G_of_u .= 0

    # Extract parameters from inputs
    S = space(component(u,1))
    f = frequency(S)
    dim = div(dimension(space(u)), dimension(S))

    S_pad = Fourier(div(N_fft,2), f)    

    u_pad = project(u, S_pad^dim)

    # Evaluate the fourier series at grid points using the IFFT
    grid_points = zeros(ComplexF64, dim, N_fft)
    for i ∈ 1:dim
        grid_points[i,:] = N_fft * FFTW.ifft(FFTW.ifftshift(component(u_pad,i).coefficients[begin:end-1]))
    end

    # Evaluate the function at the grid points
    eval_points = zeros(ComplexF64, dim, dim, N_fft)

    for j ∈ 1:N_fft
        eval_points[:,:,j] = G(grid_points[:,j])
    end

    # Use the FFT to obtain the coefficients of G(u)
    fourier_coeffs = zeros(ComplexF64, dim, dim, N_fft)
    for i ∈ 1:dim
        for j ∈ 1:dim
            fourier_coeffs[i,j,:] = 1/N_fft * FFTW.fftshift(FFTW.fft(eval_points[i,j,:]))
        end
    end

    for i ∈ 1:dim
        for j ∈ 1:dim
            component(G_of_u,i,j).coefficients[:] = project(Multiplication(Sequence(S_pad, [fourier_coeffs[i,j,:];0])), S,S).coefficients[:]
        end
    end     
end

function Newton(u::Sequence, F::Sequence, DF::LinearOperator, tol::Float64 = 1e-12, max_iter::Int64 = 20)
    count = 0;
    F!(F, u, N_fft)
    DF!(DF, u, N_fft)
    while norm(F) > tol && count <= max_iter
        u = u - DF \ F
        F!(F,u,N_fft)
        DF!(DF, u, N_fft)
        count = count + 1
    end
    return u
end

function collection_eval(time_data::Vector{Float64}, u::Sequence)
    num_points = length(time_data)
    space_data = zeros(3, num_points)

    for i ∈ 1:num_points
        space_data[:,i] = real(u(time_data[i]))
    end


    return space_data
end