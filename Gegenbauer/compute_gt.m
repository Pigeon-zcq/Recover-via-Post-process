function gt = compute_gt(lam, maxj)
% COMPUTE_GT  Precompute Pochhammer (rising factorial) ratios.
%
%   gt(j+1) = (lam)_j / j!  =  Gamma(lam+j) / (Gamma(lam) * j!)
%   for j = 0, 1, ..., maxj.
%
% Uses the stable one-step recurrence:
%   gt(j+1) = gt(j) * (lam + j - 1) / j
%
% These appear in the trigonometric expansion of Gegenbauer polynomials:
%   C_n^lam(cos theta) = sum_{j=0}^{n} gt(j) * gt(n-j) * cos((n-2j)*theta)
%
% Input:
%   lam   - Gegenbauer parameter lambda  (positive real)
%   maxj  - maximum index needed
%
% Output:
%   gt  - row vector of length maxj+1, gt(j+1) stores the j-th coefficient

    gt    = zeros(1, maxj + 1);
    gt(1) = 1.0;            % j=0: (lam)_0 / 0! = 1
    for j = 1:maxj
        gt(j+1) = gt(j) * (lam + j - 1) / j;
    end
end
