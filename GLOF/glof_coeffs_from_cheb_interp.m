function g = glof_coeffs_from_cheb_interp(xnodes, fnodes, wbar, ...
                                          m, alpha, beta, lambda, Q)
%GLOF_COEFFS_FROM_CHEB_INTERP  GLOF coefficients from Chebyshev values (Variant A).
%
%   g = GLOF_COEFFS_FROM_CHEB_INTERP(xnodes, fnodes, wbar, m, alpha, beta,
%                                    lambda, Q)
%   computes approximate GLOF coefficients
%       g_n = (1/gamma_n) * integral_0^1 (I_N f)(x) S_n(x)
%                          (-log x)^alpha x^lambda dx,
%   where I_N f is the barycentric Chebyshev interpolant defined by the
%   nodal values fnodes at the ascending Chebyshev--Lobatto nodes xnodes
%   with barycentric weights wbar.
%
%   The integral is evaluated by the same Chen--Shen Gauss--GLOF rule used
%   in GLOF_COEFFS_FROM_FOURIER, except that the Fourier partial sum is
%   replaced by I_N f.  No adaptive integral() call is used here.
%
%   Input
%       xnodes : 1-by-(N+1) Chebyshev--Lobatto nodes (ascending, in [0,1])
%       fnodes : nodal values f(xnodes)
%       wbar   : barycentric weights for xnodes
%       m      : largest GLOF degree
%       alpha  : Laguerre parameter, alpha > -1
%       beta   : GLOF parameter, beta > -1
%       lambda : weight parameter
%       Q      : number of Gauss--GLOF quadrature points
%
%   Output
%       g : (m+1)-by-1 real vector of GLOF coefficients.

    [xq, Wg, S, ~, ~, ~, gamma_vec] = ...
        glof_quadrature_matrices(Q, m, alpha, beta, lambda);

    Fq = bary_eval(xq, xnodes, fnodes, wbar);
    g = (S * (Wg(:) .* Fq(:))) ./ gamma_vec(:);
    g = real(g);
end
