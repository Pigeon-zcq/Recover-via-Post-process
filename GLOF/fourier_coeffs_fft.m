function [ell, c] = fourier_coeffs_fft(N, M, r, k)
%FOURIER_COEFFS_FFT  Discrete Fourier coefficients of f_logsing.
%
%   [ell, c] = FOURIER_COEFFS_FFT(N, M, r, k) samples
%       x_j = j/M,  j = 0,1,...,M-1,
%   computes f(x_j), takes an FFT, and extracts the (2N+1) coefficients
%   for modes ell = -N,...,N.  Choose M >> 2N+1 to approximate the
%   continuous Fourier coefficients.
%
%   Input
%       N : largest absolute mode index (>= 0)
%       M : number of FFT samples (M >= 2N+1)
%       r : exponent r > 0
%       k : non-negative integer k
%
%   Output
%       ell : (2N+1)-by-1 column vector of mode indices, -N..N
%       c   : (2N+1)-by-1 complex vector of Fourier coefficients

    xs   = (0:M-1) / M;
    fxs  = f_logsing(xs, r, k);
    chat = fft(fxs) / M;     % chat(l+1) ~ coefficient of exp(2*pi*i*l*x).

    ell = (-N:N).';
    c   = zeros(numel(ell), 1);
    for j = 1:numel(ell)
        l = ell(j);
        if l >= 0
            c(j) = chat(l + 1);
        else
            c(j) = chat(M + l + 1);
        end
    end
end
