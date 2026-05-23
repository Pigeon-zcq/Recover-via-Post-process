function hg = gcoe(hf, N, m, lm, bet, a, b, gt, nt)
% GCOE  Compute approximate Gegenbauer coefficients from Fourier coefficients.
%
% Approximates (Chen & Shu 2015, eq. 3.1):
%
%   hg(l) = (1/h_l^lm) * int_{-1}^{1} (1-y^2)^{lm-1/2} f_N(x(y)) C_l^lm(y) dy
%
% using Chebyshev-Gauss quadrature with nt-1 interior points:
%
%   int_{-1}^1 g(y)/sqrt(1-y^2) dy  ≈  (pi/nt) * sum_{i=1}^{nt-1} g(cos(pi*i/nt))
%
% so the integral becomes:
%
%   hg(l) ≈ (pi/nt)/h_l^lm * sum_{i=1}^{nt-1} C_l^lm(xi_i)*f_N(xr_i)*(1-xi_i^2)^lm
%
% where xi_i = cos(pi*i/nt) are the quadrature nodes and xr_i = transf_y2x(xi_i).
%
% Normalization constant (paper eq. 2.6):
%   h_l^lm = pi * 2^{1-2*lm} / Gamma(lm)^2 * Gamma(l+2*lm) / ((l+lm) * l!)
% Computed in log-space to avoid overflow for large lm.
%
% Inputs:
%   hf   - Fourier coefficients, hf(k+N+1) = hat_f_k, k = -N..N  [complex, length 2N+1]
%   N    - Fourier truncation order (2N+1 modes total)
%   m    - highest Gegenbauer degree to compute (l = 0,...,m)
%   lm   - Gegenbauer parameter lambda
%   bet  - coordinate transform exponent = 1/q
%   a,b  - domain endpoints
%   gt   - precomputed gt(j+1) = (lm)_j/j!  (length >= max(m,nt)+1)
%   nt   - number of Chebyshev-Gauss quadrature subintervals
%
% Output:
%   hg   - row vector of length m+1; hg(l+1) = approximate coefficient for C_l^lm

    % Chebyshev-Gauss interior nodes: xi_i = cos(pi*i/nt), i=1,...,nt-1
    xi = cos(pi * (1:nt-1) / nt);      % row vector, length nt-1

    % Map quadrature nodes to x-space
    xr = transf_y2x(xi, bet, a, b);

    % Evaluate truncated Fourier series at quadrature nodes:
    %   f_N(xr) = Re( sum_{k=-N}^{N} hat_f_k * exp(ik*pi*xr) )
    ff = zeros(1, nt-1);
    for k = -N:N
        ff = ff + real(hf(k+N+1) * exp(1i * k * pi * xr));
    end

    % log(hlm) where hlm = pi * 2^{1-2*lm} / Gamma(lm)^2
    log_hlm = log(pi) - 2*gammaln(lm) - (2*lm - 1)*log(2);

    hg = zeros(1, m + 1);
    for l = 0:m
        % C_l^lm at the nt interior quadrature nodes (gegen returns nt+1 points)
        [~, cnl_all] = gegen(l, nt, gt);
        cnl = cnl_all(2:nt);            % interior points i=1,...,nt-1, row vector

        % Quadrature sum: (pi/nt) * sum C_l^lm(xi) * f_N(xr) * (1-xi^2)^lm
        fi = cnl .* ff .* (1 - xi).^lm .* (1 + xi).^lm;

        % log(h_l^lm) = log(hlm) + log Gamma(l+2*lm) - log(l+lm) - log Gamma(l+1)
        log_hll = log_hlm + gammaln(l + 2*lm) - log(l + lm) - gammaln(l + 1);
        hg(l+1) = sum(fi) * (pi / nt) * exp(-log_hll);
    end
end
