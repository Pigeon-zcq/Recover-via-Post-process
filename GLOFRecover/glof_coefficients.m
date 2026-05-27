function c_tilde = glof_coefficients(f_j, Phi, chi_j, gamma_l)
% c_tilde = glof_coefficients(f_j, Phi, chi_j, gamma_l)
%
% Compute approximate GLOF expansion coefficients (Section 7.2 / Step 8):
%
%   c_tilde(l) = (1/gamma_l(l)) * sum_j  f_j(j) * Phi(j,l) * chi_j(j)
%
% In matrix form:  c_tilde = G^{-1} * Phi^T * W * f_j
% where  W = diag(chi_j),  G = diag(gamma_l).
%
% Inputs:
%   f_j     - (Q+1) x 1 function values at quadrature nodes
%   Phi     - (Q+1) x (M+1) GLOF basis matrix
%   chi_j   - (Q+1) x 1 GLOF quadrature weights
%   gamma_l - (M+1) x 1 normalization constants
%
% Output:
%   c_tilde - (M+1) x 1 GLOF coefficients

weighted_f = chi_j(:) .* f_j(:);        % W * f_j
c_tilde    = (Phi' * weighted_f) ./ gamma_l(:);
end
