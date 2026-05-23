function ermax = run_reconstruction(N, lam, m, p, q, a, b, f_func, f_exact_func)
% RUN_RECONSTRUCTION  Single-call Gegenbauer reconstruction driver.
%
% Computes Fourier coefficients of f_func, performs Gegenbauer reconstruction
% with parameters (N, lambda=lam, m) on domain [a,b], and returns the
% L_inf error against f_exact_func at N+1 Chebyshev output nodes.
%
% Inputs:
%   N           - Fourier truncation order
%   lam         - Gegenbauer parameter lambda
%   m           - number of Gegenbauer terms (l = 0,...,m)
%   p, q        - singularity exponent s = p/q
%   a, b        - domain endpoints
%   f_func      - function handle for f (used to compute Fourier coefficients)
%   f_exact_func- function handle for the exact solution (for error)
%
% Output:
%   ermax  - L_inf error at N+1 Chebyshev nodes

    bet = 1/q;
    nt  = 3*N;
    nzn = N;

    hf      = getPointValuefromFC(N, f_func);
    gt      = compute_gt(lam, max(m, nt) + 1);
    hg      = gcoe(hf, N, m, lam, bet, a, b, gt, nt);
    [u, xx] = gegenbauer_reconstruct(hg, m, lam, bet, a, b, nzn, gt);

    ex    = arrayfun(f_exact_func, xx);
    ermax = max(abs(u - ex));
end
