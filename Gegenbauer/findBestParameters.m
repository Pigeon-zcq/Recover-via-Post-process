function [mbest, lbest, err_best] = findBestParameters(N, p, q, a, b, f_func, f_exact_func)
% FINDBESTPARAMETERS  Grid search for optimal Gegenbauer parameters (lambda, m).
%
% Searches all pairs (lambda, m) with lambda in [1, N/2] and m in [0, N/2],
% and returns the pair achieving the smallest L_inf reconstruction error.
%
% For each lambda:
%   - Fourier coefficients are already computed once (shared across all m).
%   - Gegenbauer coefficients hg(0),...,hg(N/2) are computed once.
%   - The reconstruction sum is built incrementally as m increases, so each
%     (lambda, m) pair is evaluated at O(N) cost rather than O(N^2).
%
% Inputs:
%   N           - Fourier truncation order
%   p, q        - singularity exponent s = p/q
%   a, b        - domain endpoints
%   f_func      - function handle for f  (for Fourier coefficient computation)
%   f_exact_func- function handle for exact f  (for error computation)
%
% Outputs:
%   mbest    - optimal m
%   lbest    - optimal lambda
%   err_best - corresponding L_inf error

    bet  = 1/q;
    nt   = 3*N;
    nzn  = N;
    mmax = floor(N/2);

    hf = getPointValuefromFC(N, f_func);

    mbest    = 0;
    lbest    = 1;
    err_best = Inf;
    ex_at_nodes = [];   % set on first iteration

    for lam = 1 : floor(N/2)
        gt = compute_gt(lam, max(mmax, nt) + 1);

        % Evaluate f_N at the nt quadrature nodes (same for all m with this lam)
        xi = cos(pi * (1:nt-1) / nt);
        xr = transf_y2x(xi, bet, a, b);
        ff = zeros(1, nt-1);
        for k = -N:N
            ff = ff + real(hf(k+N+1) * exp(1i*k*pi*xr));
        end
        log_hlm = log(pi) - 2*gammaln(lam) - (2*lam - 1)*log(2);

        % Accumulate Gegenbauer sum incrementally over m = 0, 1, ..., mmax
        u_accum = zeros(1, nzn+1);

        for l = 0:mmax
            % Coefficient hg(l) for the current lambda
            [~, cnl_q] = gegen(l, nt, gt);
            cnl_int    = cnl_q(2:nt);
            fi         = cnl_int .* ff .* (1 - xi).^lam .* (1 + xi).^lam;
            log_hll    = log_hlm + gammaln(l+2*lam) - log(l+lam) - gammaln(l+1);
            hg_l       = sum(fi) * (pi / nt) * exp(-log_hll);

            % Gegenbauer polynomial at output nodes
            [yy_l, cnl_r] = gegen(l, nzn, gt);
            if lam == 1 && l == 0
                % Output nodes are fixed (depend only on nzn)
                xx_nodes    = transf_y2x(yy_l, bet, a, b);
                ex_at_nodes = arrayfun(f_exact_func, xx_nodes);
            end
            u_accum = u_accum + hg_l * cnl_r;

            err_now = max(abs(u_accum - ex_at_nodes));
            if err_now < err_best
                err_best = err_now;
                lbest    = lam;
                mbest    = l;
            end
        end
    end
end
