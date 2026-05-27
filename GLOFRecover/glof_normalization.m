function gamma_l = glof_normalization(alpha, beta, M)
% gamma_l = glof_normalization(alpha, beta, M)
%
% Compute GLOF orthogonality normalization constants (Section 3.3 / Step 7):
%
%   gamma_l^(alpha,beta) = Gamma(l+alpha+1) / [ (beta+1)^(alpha+1) * Gamma(l+1) ]
%
% Uses log-gamma to avoid overflow for large l or alpha.
%
% Inputs:
%   alpha - Laguerre weight parameter
%   beta  - GLOF log-mapping parameter
%   M     - highest GLOF order
%
% Output:
%   gamma_l - (M+1) x 1 vector

l = (0 : M)';
log_gamma_l = gammaln(l + alpha + 1) ...
            - gammaln(l + 1)         ...
            - (alpha + 1) * log(beta + 1);
gamma_l = exp(log_gamma_l);
end
