function f_j = fourier_partial_sum(f_hat, x_j, NF)
% f_j = fourier_partial_sum(f_hat, x_j, NF)
%
% Evaluate the Fourier partial sum at non-uniform nodes x_j (Section 5):
%
%   f_j(i) = sum_{n=-NF}^{NF} f_hat(n) * exp(2*pi*i*n*x_j(i))
%
% f_hat is ordered as n = -NF, -NF+1, ..., NF  (length 2*NF+1).
% Since x_j are non-uniform, direct matrix multiplication is used (O(Q*NF)).
%
% Inputs:
%   f_hat - (2*NF+1) x 1 complex Fourier coefficients
%   x_j   - (Q+1) x 1 quadrature nodes in (0,1)
%   NF    - Fourier half-bandwidth
%
% Output:
%   f_j - (Q+1) x 1 complex vector (take real() for real-valued f)

x_j   = x_j(:);
n_vec = -NF : NF;                           % 1 x (2*NF+1)
E     = exp(2i * pi * x_j * n_vec);         % (Q+1) x (2*NF+1)
f_j   = E * f_hat(:);
end
