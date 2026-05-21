function x = cheb_lobatto_nodes(N)
%CHEB_LOBATTO_NODES  Chebyshev--Lobatto nodes on [0,1] in ascending order.
%
%   x = CHEB_LOBATTO_NODES(N) returns a 1-by-(N+1) row vector containing
%       x_j = (1 + cos(pi*j/N))/2,    j = 0,1,...,N,
%   re-ordered so that x is monotonically increasing from 0 to 1.
%
%   Input
%       N : polynomial degree (N >= 1)
%
%   Output
%       x : 1-by-(N+1) ascending vector with x(1) = 0, x(end) = 1.

    j = 0:N;
    x = (1 + cos(pi*j/N)) / 2;
    x = fliplr(x);
    x(1)   = 0;
    x(end) = 1;
end
