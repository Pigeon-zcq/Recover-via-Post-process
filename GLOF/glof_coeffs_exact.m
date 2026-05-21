function g = glof_coeffs_exact(f_handle, m, alpha, beta, lambda, Q)
%GLOF_COEFFS_EXACT  GLOF coefficients of a function from exact values.
%
%   g = GLOF_COEFFS_EXACT(f_handle, m, alpha, beta, lambda, Q) returns
%   the GLOF coefficients of f computed by the same Chen--Shen Gauss--GLOF
%   quadrature used for post-processing, applied directly to f (NOT to a
%   truncated spectral expansion).
%
%   This is the generic direct-projection routine for an arbitrary function
%   handle.  For high-order projections of the analytic model
%
%       f(x) = x^r * (-log x)^k,
%
%   prefer GLOF_COEFFS_LOGSING_STABLE(r,k,m,alpha,beta,lambda,Q).  That
%   specialized routine evaluates the transformed coefficient integral in a
%   more stable way and avoids catastrophic cancellation in high-order
%   Laguerre quadrature.  This routine does not try to infer r and k from
%   f_handle, because a general function handle cannot be reliably parsed.
%
%   Input
%       f_handle : function handle accepting an array x of values in [0,1]
%                  and returning f(x).
%       m        : largest GLOF degree
%       alpha    : Laguerre parameter
%       beta     : GLOF parameter
%       lambda   : weight parameter
%       Q        : number of Gauss--GLOF quadrature points
%
%   Output
%       g : (m+1)-by-1 real vector of exact GLOF coefficients.

    [xq, Wg, S, ~, ~, ~, gamma_vec] = ...
        glof_quadrature_matrices(Q, m, alpha, beta, lambda);

    Fq = f_handle(xq);
    g = (S * (Wg(:) .* Fq(:))) ./ gamma_vec(:);
    g = real(g);
end
