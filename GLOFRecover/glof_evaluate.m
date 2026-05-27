function R = glof_evaluate(c_tilde, z_eval, alpha, beta, lambda)
% R = glof_evaluate(c_tilde, z_eval, alpha, beta, lambda)
%
% Evaluate the GLOF reconstruction at arbitrary points z_eval (Section 8):
%
%   R(i) = sum_{l=0}^{M} c_tilde(l) * z_i^nu * L_l^(alpha)(Y_i)
%
%   where  Y_i = -(beta+1)*log(z_i),   nu = (beta-lambda)/2.
%
% In matrix form:  R = Psi * c_tilde,   Psi(i,l+1) = z_i^nu * L_l^(alpha)(Y_i).
%
% Inputs:
%   c_tilde - (M+1) x 1 GLOF coefficients
%   z_eval  - N_eval x 1 evaluation points in (0,1)
%   alpha   - Laguerre parameter
%   beta    - GLOF log-mapping parameter
%   lambda  - GLOF weight exponent
%
% Output:
%   R - N_eval x 1 reconstruction values

nu     = (beta - lambda) / 2;
z_eval = z_eval(:);
Y_eval = -(beta + 1) * log(z_eval);         % mapped variable
z_nu   = exp(nu * log(z_eval));              % z^nu, safe for small z

M     = length(c_tilde) - 1;
L_mat = laguerre_poly_matrix(Y_eval, alpha, M);   % N_eval x (M+1)
Psi   = z_nu .* L_mat;                            % evaluation matrix

R = Psi * c_tilde(:);
end
