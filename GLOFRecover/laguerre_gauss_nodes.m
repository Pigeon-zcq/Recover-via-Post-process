function [y, w] = laguerre_gauss_nodes(Q, alpha)
% [y, w] = laguerre_gauss_nodes(Q, alpha)
%
% Compute Q+1 generalized Laguerre-Gauss quadrature nodes y and weights w
% via the Golub-Welsch algorithm.
%
% The rule integrates x^alpha * exp(-x) * p(x) exactly for p in P_{2Q+1}.
%
% Inputs:
%   Q     - highest node index; Q+1 nodes total
%   alpha - Laguerre weight parameter, alpha > -1
%
% Outputs:
%   y  - (Q+1) x 1 vector of nodes (positive)
%   w  - (Q+1) x 1 vector of weights (positive)

n = Q + 1;   % number of nodes

% Diagonal entries of Jacobi matrix: a_k = 2k + alpha + 1, k = 0,...,n-1
k  = (0 : n-1)';
a  = 2*k + alpha + 1;

% Off-diagonal entries: sqrt(beta_k) = sqrt(k*(k+alpha)), k = 1,...,n-1
k2 = (1 : n-1)';
b  = sqrt(k2 .* (k2 + alpha));

% Symmetric tridiagonal Jacobi matrix
T = diag(a) + diag(b, 1) + diag(b, -1);

% Eigendecomposition (MATLAB uses stable symmetric solver)
[V, D] = eig(T);
[y, idx] = sort(diag(D));
V = V(:, idx);

% Golub-Welsch weights: w_j = Gamma(alpha+1) * V(1,j)^2
w = gamma(alpha + 1) * V(1, :)' .^ 2;
end
