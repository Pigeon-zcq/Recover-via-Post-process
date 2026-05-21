function w = bary_weights_cheb_lobatto(N)
%BARY_WEIGHTS_CHEB_LOBATTO  Barycentric weights for Chebyshev--Lobatto nodes.
%
%   w = BARY_WEIGHTS_CHEB_LOBATTO(N) returns the (1-by-(N+1)) row vector
%   of barycentric weights associated with the Chebyshev--Lobatto nodes
%   returned by CHEB_LOBATTO_NODES(N).  The convention used here matches
%   the ascending ordering of the nodes (from 0 to 1).
%
%   Input
%       N : polynomial degree (N >= 1)
%
%   Output
%       w : 1-by-(N+1) row vector of barycentric weights.

    w      = ones(1, N+1);
    w(1)   = 0.5;
    w(end) = 0.5;
    w      = w .* (-1).^(0:N);
    w      = fliplr(w);
end
