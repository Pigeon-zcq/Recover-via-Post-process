function g = glof_coeffs_from_cheb_weighted(xnodes, fnodes, wbar, ...
                                            m, alpha, beta, lambda, Q)
%GLOF_COEFFS_FROM_CHEB_WEIGHTED  Compatibility wrapper for Chebyshev data.
%
%   This entry point is kept for older scripts that used the former
%   "Variant B" name.  The current experiment specification computes GLOF
%   coefficients from Chebyshev collocation values by first forming the
%   barycentric interpolant F_N = I_N f and then applying the Chen--Shen
%   Gauss--GLOF quadrature rule.  Therefore this wrapper delegates to
%   GLOF_COEFFS_FROM_CHEB_INTERP and does not call MATLAB integral().
%
%   Q is optional for backward compatibility.  When omitted, a conservative
%   quadrature order based on the Chebyshev degree is used.

    if nargin < 8 || isempty(Q)
        N = numel(xnodes) - 1;
        Q = max(60, 3*N);
    end

    g = glof_coeffs_from_cheb_interp(xnodes, fnodes, wbar, ...
                                     m, alpha, beta, lambda, Q);
end
