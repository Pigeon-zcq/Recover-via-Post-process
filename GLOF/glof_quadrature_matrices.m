function [xq, Wg, S, L, y, wLag, gamma_vec] = glof_quadrature_matrices(Q, m, alpha, beta, lambda)
%GLOF_QUADRATURE_MATRICES  Chen--Shen Gauss--GLOF quadrature data.
%
%   [xq, Wg, S, L, y, wLag, gamma_vec] =
%       GLOF_QUADRATURE_MATRICES(Q, m, alpha, beta, lambda)
%   returns Q mapped quadrature nodes and weights for
%
%       integral_0^1 h(x) (-log x)^alpha x^lambda dx
%           ~= sum_q h(xq(q)) Wg(q).
%
%   The construction follows the log mapping y = -(beta+1) log(x):
%
%       xq = exp(-yq/(beta+1)),
%       Wg = xq^(lambda-beta) (beta+1)^(-alpha-1) wLag,
%
%   where yq,wLag are generalized Gauss--Laguerre nodes and weights for
%   y^alpha exp(-y).  S contains the GLOF basis values
%   S_n(xq) = xq^theta L_n^{(alpha)}(yq), theta=(beta-lambda)/2.
%
%   This is the quadrature rule used for GLOF coefficient post-processing
%   from Fourier coefficients or Chebyshev collocation values.

    theta = (beta - lambda) / 2;

    [y, wLag] = gauss_laguerre(Q, alpha);
    xq = exp(-y / (beta + 1));

    Wg = (xq.^(lambda - beta)) ...
       .* (beta + 1)^(-alpha - 1) ...
       .* wLag;

    L = laguerre_matrix(m, alpha, y);
    S = bsxfun(@times, xq(:).'.^theta, L);

    gamma_vec = glof_norm_const((0:m).', alpha, beta);
end
