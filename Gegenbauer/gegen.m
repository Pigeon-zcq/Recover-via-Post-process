function [yy, cnl] = gegen(n, nn, gt)
% GEGEN  Evaluate Gegenbauer polynomial C_n^lambda at nn+1 Chebyshev nodes.
%
% Computes C_n^lambda(cos(theta_i)) for theta_i = pi*i/nn, i = 0,...,nn,
% using the trigonometric series representation:
%
%   C_n^lambda(cos theta) = sum_{j=0}^{n} gt(j)*gt(n-j)*cos((n-2j)*theta)
%
% where gt(j) = (lambda)_j / j!  (precomputed via compute_gt).
%
% Inputs:
%   n   - polynomial degree
%   nn  - number of subintervals; generates nn+1 Chebyshev nodes
%   gt  - precomputed Pochhammer ratios (1-indexed: gt(j+1) = (lam)_j/j!)
%         must have length >= n+1
%
% Outputs:
%   yy (1:nn+1)  = cos(pi*i/nn)         Chebyshev nodes in y-space
%   cnl(1:nn+1)  = C_n^lambda(yy)       polynomial values at those nodes

    yy  = zeros(1, nn + 1);
    cnl = zeros(1, nn + 1);
    for i = 0:nn
        theta    = pi * i / nn;
        yy(i+1)  = cos(theta);
        s = 0;
        for j = 0:n
            s = s + gt(j+1) * gt(n-j+1) * cos(theta * (n - 2*j));
        end
        cnl(i+1) = s;
    end
end
