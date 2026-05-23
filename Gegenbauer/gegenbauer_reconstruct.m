function [u, xx] = gegenbauer_reconstruct(hg, m, lam, bet, a, b, nzn, gt)
% GEGENBAUER_RECONSTRUCT  Evaluate Gegenbauer approximation at Chebyshev nodes.
%
% Computes:
%   u(x_i) = sum_{l=0}^{m} hg(l) * C_l^lambda(y(x_i))
%
% at nzn+1 output nodes x_i = transf_y2x(cos(pi*i/nzn)), i=0,...,nzn.
%
% Inputs:
%   hg   - Gegenbauer coefficients, hg(l+1) for l=0,...,m  (from gcoe)
%   m    - highest degree
%   lam  - Gegenbauer parameter lambda
%   bet  - coordinate transform exponent = 1/q
%   a,b  - domain endpoints
%   nzn  - number of output subintervals (nzn+1 output points)
%   gt   - precomputed gt(j+1) = (lam)_j/j!  (length >= m+1)
%
% Outputs:
%   u    - reconstructed values, length nzn+1
%   xx   - corresponding x-coordinates, length nzn+1

    u  = zeros(1, nzn + 1);
    xx = zeros(1, nzn + 1);

    for l = 0:m
        [yy_l, cnl_l] = gegen(l, nzn, gt);
        if l == 0
            % All gegen calls share the same yy (depends only on nzn)
            xx = transf_y2x(yy_l, bet, a, b);
        end
        u = u + hg(l+1) * cnl_l;
    end
end
