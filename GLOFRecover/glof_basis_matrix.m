function Phi = glof_basis_matrix(x_j, y_j, alpha, beta, lambda, M)
% Phi = glof_basis_matrix(x_j, y_j, alpha, beta, lambda, M)
%
% Compute the GLOF basis matrix at quadrature nodes (Section 6):
%
%   Phi(j, l+1) = x_j^nu * L_l^(alpha)(y_j),   nu = (beta - lambda) / 2
%
% x_j^nu is evaluated as exp(nu * log(x_j)) for numerical safety when x_j is small.
%
% Inputs:
%   x_j    - (Q+1) x 1 quadrature nodes in (0,1)
%   y_j    - (Q+1) x 1 Laguerre-Gauss nodes (same index set)
%   alpha  - Laguerre parameter
%   beta   - GLOF log-mapping parameter
%   lambda - GLOF weight exponent
%   M      - highest GLOF order
%
% Output:
%   Phi - (Q+1) x (M+1) matrix

nu    = (beta - lambda) / 2;
x_j   = x_j(:);
x_nu  = exp(nu * log(x_j));             % x_j^nu, safe for small x_j

L_mat = laguerre_poly_matrix(y_j, alpha, M);   % (Q+1) x (M+1)
Phi   = x_nu .* L_mat;                        % broadcast column x_nu
end
