function [z, w] = gauss_laguerre(Q, alpha)
%GAUSS_LAGUERRE  Generalized Gauss--Laguerre nodes and weights.
%
%   [z, w] = GAUSS_LAGUERRE(Q, alpha) returns Q nodes z and weights w
%   such that, for sufficiently regular h,
%       integral_0^inf h(z) * z^alpha * exp(-z) dz  ~=  sum_q w(q)*h(z(q)).
%
%   The Jacobi (tridiagonal) matrix associated with the three-term
%   recurrence of the generalized Laguerre polynomials is built and its
%   eigendecomposition is used (Golub--Welsch algorithm).
%
%   Input
%       Q     : number of quadrature points (Q >= 1)
%       alpha : Laguerre parameter, alpha > -1
%
%   Output
%       z : Q-by-1 column vector of nodes in (0, inf)
%       w : Q-by-1 column vector of weights

    j = (0:Q-1).';
    diag_main = 2*j + alpha + 1;
    diag_off  = sqrt((1:Q-1).' .* ((1:Q-1).' + alpha));

    J = diag(diag_main) + diag(diag_off, 1) + diag(diag_off, -1);

    [V, D] = eig(J);
    z = diag(D);
    [z, idx] = sort(z);
    V = V(:, idx);

    % Weights: gamma(alpha+1) * (first eigenvector component)^2.
    w = gamma(alpha + 1) * (V(1, :).').^2;
end
