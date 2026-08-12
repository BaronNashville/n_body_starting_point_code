function G_of_u_pt2pt!(G_of_u::Sequence, u::Sequence, G::Function; N_fft = 2^14)
    # INPUTS:
    # 1) u  = {u_n}_{n=-N}^{N} sequence of fourier coefficients
    # 2) g : R → R analytic function

    # Set ouput to 0
    G_of_u.coefficients .= 0

    # Extract parameters from inputs
    f = frequency(space(G_of_u))
    N_out = order(G_of_u)

    S_pad = Fourier(div(N_fft, 2), f)
    u_pad = project(u, S_pad)

    # Evaluate the fourier series at grid points using the IFFT
    grid_points = Sequence(S_pad, [FFTW.fft(FFTW.ifftshift(u_pad.coefficients[begin:end-1])); 0])

    # Evaluate the function at the grid points
    eval_points = zeros(ComplexF64, S_pad)

    for n ∈ -div(N_fft, 2):div(N_fft, 2)-1
        eval_points[n] = G(grid_points[n])
    end

    # Use the FFT to obtain the coefficients of G(u)
    G_of_u.coefficients[:] = project(Sequence(S_pad, [FFTW.fftshift(FFTW.ifft(eval_points.coefficients[begin:end-1])); 0]), Fourier(N_out, f)).coefficients[:]
end

function G_of_u_vec2pt!(G_of_u::Sequence, u::Sequence, G::Function; N_fft = 2^14)
    # INPUTS:
    # 1) u  = {u_n}_{n=-N}^{N} sequence of fourier coefficients with n components
    # 2) g : Rⁿ → R analytic function

    # Set ouput to 0
    G_of_u.coefficients .= 0

    # Extract parameters from inputs
    f = frequency(space(G_of_u))[1]
    N_out = order(G_of_u)
    dim = size(N_out)

    S_pad = Fourier(div(N_fft, 2), f)
    u_pad = project(u, S_pad^dim)

    # Evaluate the fourier series at grid points using the IFFT
    grid_points = zeros(ComplexF64, dim, N_fft)
    for i ∈ 1:dim
        grid_points[i, :] = FFTW.fft(FFTW.ifftshift(block(u_pad, i).coefficients[begin:end-1]))
    end

    # Evaluate the function at the grid points
    eval_points = zeros(ComplexF64, N_fft)

    for j ∈ 1:N_fft
        eval_points[j] = G(grid_points[:, j])
    end

    # Use the FFT to obtain the coefficients of G(u)
    fourier_coeffs = FFTW.fftshift(FFTW.ifft(eval_points))

    G_of_u.coefficients[:] = project(Sequence(S_pad, [fourier_coeffs; 0]), Fourier(N_out[1], f)).coefficients[:]
end

function G_of_u_vec2vec!(G_of_u::Sequence, u::Sequence, G::Function; N_fft = 2^14)
    # INPUTS:
    # 1) u  = {u_n}_{n=-N}^{N} sequence of fourier coefficients with n components
    # 2) g : Rⁿ → Rⁿ analytic function

    # Set ouput to 0
    G_of_u.coefficients .= 0

    # Extract parameters from inputs
    f = frequency(space(G_of_u))[1]
    N_out = order(G_of_u)
    dim = size(N_out)

    S_pad = Fourier(div(N_fft, 2), f)
    u_pad = project(u, S_pad^dim)

    # Evaluate the fourier series at grid points using the IFFT
    grid_points = zeros(ComplexF64, dim, N_fft)
    for i ∈ 1:dim
        grid_points[i, :] = FFTW.fft(FFTW.ifftshift(block(u_pad, i).coefficients[begin:end-1]))
    end

    # Evaluate the function at the grid points
    eval_points = zeros(ComplexF64, dim, N_fft)

    for j ∈ 1:N_fft
        eval_points[:, j] = G(grid_points[:, j])
    end

    # Use the FFT to obtain the coefficients of G(u)
    fourier_coeffs = zeros(ComplexF64, dim, N_fft)
    for i ∈ 1:dim
        fourier_coeffs[i, :] = FFTW.fftshift(FFTW.ifft(eval_points[i, :]))
    end

    for i ∈ 1:dim
        block(G_of_u, i).coefficients[:] = project(Sequence(S_pad, [fourier_coeffs[i, :]; 0]), Fourier(N_out[i], f)).coefficients[:]
    end
end

function G_of_u_vec2mat!(G_of_u::LinearOperator, u::Sequence, G::Function; N_fft = 2^14)
    # INPUTS:
    # 1) u  = {u_n}_{n=-N}^{N} sequence of fourier coefficients with n components
    # 2) G : Rⁿ → Rⁿ x Rⁿ  analytic function

    # Set ouput to 0
    G_of_u.coefficients .= 0

    # Extract parameters from inputs
    f = frequency(space(G_of_u))[1][1]
    N_out = order(G_of_u)
    dim = size(N_out)

    S_pad = Fourier(div(N_fft, 2), f)
    u_pad = project(u, S_pad^dim)

    # Evaluate the fourier series at grid points using the IFFT
    grid_points = zeros(ComplexF64, dim, N_fft)
    for i ∈ 1:dim
        grid_points[i, :] = FFTW.fft(FFTW.ifftshift(block(u_pad, i).coefficients[begin:end-1]))
    end

    # Evaluate the function at the grid points
    eval_points = zeros(ComplexF64, dim, dim, N_fft)

    for j ∈ 1:N_fft
        eval_points[:, :, j] = G(grid_points[:, j])
    end

    # Use the FFT to obtain the coefficients of G(u)
    fourier_coeffs = zeros(ComplexF64, dim, dim, N_fft)
    for i ∈ 1:dim
        for j ∈ 1:dim
            fourier_coeffs[i, j, :] = FFTW.fftshift(FFTW.ifft(eval_points[i, j, :]))
        end
    end

    for i ∈ 1:dim
        for j ∈ 1:dim
            block(G_of_u, i, j).coefficients[:] = project(Multiplication(Sequence(S_pad, [fourier_coeffs[i, j, :]; 0])), Fourier(N_out[j], f), Fourier(N_out[j], f)).coefficients[:]
        end
    end
end

function cheb_G_of_u_pt2pt!(G_of_u::Sequence, u::Sequence, G::Function;  N_fft = 2^14, N_out = order(G_of_u))
    # INPUTS:
    # 1) u = {u_n}_{n = 0}^{N} sequence of Chebyshev coefficients
    # 2) G : R → R analytic function

    N_in = order(u)

    # Set output to 0
    G_of_u.coefficients .= 0

    # Step 1: Use the natural transformation to go form Chebyshev
    # series to Fourier series t → cosθ
    fourier_coeffs = [reverse(u.coefficients[begin+1:end], u.coefficients[:])]
    fourier_sequence = Sequence(Fourier(N_in,1.0), fourier_coeffs)

    # Step 2: Use the methods above to compute the Fourier coefficients of G(u(cosθ))
    transformed_fourier_sequence = zeros(Sequence(Fourier(N_out, 1.0)))
    G_of_u!(transformed_fourier_sequence, fourier_sequence, G)

    # Step 3: Return to Chebyshev world
    G_of_u.coefficients[:] = transformed_fourier_sequence[0:end]
end
