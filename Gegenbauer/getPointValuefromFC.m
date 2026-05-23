function hf = getPointValuefromFC(N, f_func)
% GETPOINTVALUEFROMFC  Compute Fourier coefficients by numerical integration.
%
% Approximates:
%   hat_f_k = 0.5 * int_{-1}^{1} f(x) * exp(-i*k*pi*x) dx,  k = -N,...,N
%
% using 3-point Gaussian quadrature on M = min(100*N, 4000) panels of
% width 2/M each, covering the interval [-1, 1]:
%
%   Gauss nodes per panel:  xm - x1*dx,  xm,  xm + x1*dx
%   Gauss weights:           wt1,         wt0,  wt1
%   x1 = sqrt(3/5),  wt1 = 5/9,  wt0 = 8/9   (scaled by panel half-width dx)
%
% Inputs:
%   N      - Fourier truncation order
%   f_func - function handle  f_func(x) evaluating f at a scalar x
%
% Output:
%   hf  - complex row vector of length 2N+1
%         hf(k+N+1) = hat_f_k for k = -N,...,N
%         Satisfies hf(-k) = conj(hf(k)) for real-valued f.

    x1  = sqrt(0.6);            % = sqrt(3/5), Gauss quadrature node offset
    wt1 = 5.0/9.0;              % endpoint Gauss weight
    wt0 = 8.0/9.0;              % midpoint Gauss weight
    M   = min(100*N, 4000);    % total number of panels
    dx  = 1.0 / M;             % half-panel width (full panel = 2/M)
    wt1 = wt1 * dx;            % scale weights by panel half-width
    wt0 = wt0 * dx;

    % Panel center indices: odd integers from -M+1 to M-1
    j_vals = (-M+1 : 2 : M-1);    % length M

    hf_pos = zeros(1, N+1);       % hf_pos(k+1) = hat_f_k for k = 0,...,N
    for k = 0:N
        xn = (j_vals - x1) * dx;
        xm =  j_vals        * dx;
        xp = (j_vals + x1) * dx;

        fn = arrayfun(f_func, xn) .* exp(-1i*k*pi*xn);
        fm = arrayfun(f_func, xm) .* exp(-1i*k*pi*xm);
        fp = arrayfun(f_func, xp) .* exp(-1i*k*pi*xp);

        hf_pos(k+1) = 0.5 * sum(wt1*(fn + fp) + wt0*fm);
    end

    % Fill the full array using conjugate symmetry
    hf = zeros(1, 2*N+1);
    for k = 0:N
        hf(k+N+1)  = hf_pos(k+1);
        hf(-k+N+1) = conj(hf_pos(k+1));
    end
end
