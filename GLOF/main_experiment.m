%MAIN_EXPERIMENT  Master driver for the GLOF reconstruction experiment.
%
%   Loops over several (r, k) parameter combinations and several truncation
%   sizes N.  For each (r, k, N) it records:
%       EF   : traditional Fourier maximum error
%       EC   : traditional Chebyshev interpolation maximum error
%       EGF  : GLOF reconstruction from Fourier coefficients
%       EGC  : GLOF reconstruction from Chebyshev values
%       QrelF/QrelC : ||g_Q-g_{2Q}||_2 / ||g_{2Q}||_2 diagnostics
%
%   Two GLOF parameter strategies are available:
%       mode = 'fixed'        :  alpha=0, beta=2r, lambda=0
%       mode = 'Ndependent'   :  alpha=aN, beta=2r+bN, lambda=bN
%
%   The Chen--Shen style 'Ndependent' setting is the more interesting one
%   for proving rapid recovery from standard spectral data.
%
%   At the end, a multi-panel figure summarises the convergence behaviour
%   across cases.

clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(here);

mode = 'Ndependent';   % or 'fixed'
cases = [
    struct('r', 0.5, 'k', 1)
];

Ns = [16 32 64 128 256];

xref = linspace(0, 1, 5001);

EF   = zeros(numel(cases), numel(Ns));
EC   = zeros(numel(cases), numel(Ns));
EGF  = zeros(numel(cases), numel(Ns));
EGC   = zeros(numel(cases), numel(Ns));
QrelF = zeros(numel(cases), numel(Ns));
QrelC = zeros(numel(cases), numel(Ns));

for ic = 1:numel(cases)
    r = cases(ic).r;
    k = cases(ic).k;
    fprintf('\n>>> Case %d: r = %.2f, k = %d, mode = %s\n', ic, r, k, mode);
    f_handle = @(x) f_logsing(x, r, k);
    fref = f_handle(xref);

    for iN = 1:numel(Ns)
        N = Ns(iN);

        % ----- GLOF parameters --------------------------------------------
        % For the pure model x^r (-log x)^k the tuned GLOF expansion has
        % only k+1 nonzero coefficients.  Keeping m = k avoids reconstructing
        % high GLOF modes that are driven only by Fourier/Chebyshev data noise
        % and become extremely ill-conditioned near x = 0.
        switch mode
            case 'fixed'
                alpha  = 0;
                beta   = 2*r;
                lambda = 0;
                m      = k;
                Q      = 2*N;
            case 'Ndependent'
                a_par = 0.05; b_par = 0.05;c_par = 0.01;
                alpha  = a_par * N;
                beta   = 2*r + b_par * N;
                lambda = b_par * N;
                m      = ceil(k + c_par * N);
                Q      = 2*N;
            otherwise
                error('Unknown mode "%s"', mode);
        end

        % ----- Fourier baseline --------------------------------------------
        ell = (-N:N).';
        c   = fourier_coeffs_quad(f_handle, ell);
        fF  = eval_fourier_sum(xref, ell, c);
        EF(ic, iN) = max(abs(fref - fF));

        % ----- GLOF from Fourier coefficients ------------------------------
        gF  = glof_coeffs_from_fourier(ell, c, m, alpha, beta, lambda, Q);
        gF2 = glof_coeffs_from_fourier(ell, c, m, alpha, beta, lambda, 2*Q);
        RF  = eval_glof_reconstruction(xref, gF, alpha, beta, lambda);
        EGF(ic, iN) = max(abs(fref - RF));
        QrelF(ic, iN) = norm(gF - gF2) / (norm(gF2) + 1e-300);

        % ----- Chebyshev baseline ------------------------------------------
        xn   = cheb_lobatto_nodes(N);
        fn   = f_logsing(xn, r, k);
        wbar = bary_weights_cheb_lobatto(N);
        fC   = bary_eval(xref, xn, fn, wbar);
        EC(ic, iN) = max(abs(fref - fC));

        % ----- GLOF from Chebyshev values ----------------------------------
        gC   = glof_coeffs_from_cheb_interp(xn, fn, wbar, ...
                                            m, alpha, beta, lambda, Q);
        gC2  = glof_coeffs_from_cheb_interp(xn, fn, wbar, ...
                                            m, alpha, beta, lambda, 2*Q);
        RC   = eval_glof_reconstruction(xref, gC, alpha, beta, lambda);
        EGC(ic, iN) = max(abs(fref - RC));
        QrelC(ic, iN) = norm(gC - gC2) / (norm(gC2) + 1e-300);

        fprintf(['  N=%4d | Fourier %.2e | Cheb %.2e | ' ...
                 'GLOF-F %.2e | GLOF-C %.2e | QrelF %.2e | QrelC %.2e\n'], ...
                N, EF(ic,iN), EC(ic,iN), EGF(ic,iN), EGC(ic,iN), ...
                QrelF(ic,iN), QrelC(ic,iN));
    end
end

% ----- Convergence summary plot ------------------------------------------

figure('Name', sprintf('GLOF reconstruction summary (mode=%s)', mode), ...
       'Position', [80 80 1100 700]);

for ic = 1:numel(cases)
    subplot(1, numel(cases), ic);
    semilogy(Ns, EF(ic,:),   '-o', 'LineWidth', 1.4); hold on;
    semilogy(Ns, EC(ic,:),   '-s', 'LineWidth', 1.4);
    semilogy(Ns, EGF(ic,:),  '-^', 'LineWidth', 1.4);
    semilogy(Ns, EGC(ic,:),  '-d', 'LineWidth', 1.4);
    legend('Fourier','Chebyshev','GLOF-F','GLOF-C', 'Location','best');
    grid on;
    xlabel('N'); ylabel('max error');
    title(sprintf('r=%.2f, k=%d', cases(ic).r, cases(ic).k));
end

sgtitle(sprintf('GLOF reconstruction errors (parameters: %s)', mode));

figure('Name', sprintf('Gauss-GLOF quadrature diagnostics (mode=%s)', mode), ...
       'Position', [120 120 900 500]);
for ic = 1:numel(cases)
    subplot(1, numel(cases), ic);
    semilogy(Ns, QrelF(ic,:), '-o', 'LineWidth', 1.4); hold on;
    semilogy(Ns, QrelC(ic,:), '-s', 'LineWidth', 1.4);
    grid on;
    xlabel('N'); ylabel('relative coefficient difference');
    title(sprintf('r=%.2f, k=%d', cases(ic).r, cases(ic).k));
    legend('Fourier data', 'Cheb data', 'Location','best');
end

% ----- Print final summary -----------------------------------------------

fprintf('\n=== Summary table (mode = %s) ===\n', mode);
for ic = 1:numel(cases)
    fprintf('\nCase r=%.2f, k=%d\n', cases(ic).r, cases(ic).k);
    fprintf('  N         Fourier      Cheb       GLOF-F     GLOF-C     QrelF      QrelC\n');
    for iN = 1:numel(Ns)
        fprintf('  %4d   %.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n', ...
                Ns(iN), EF(ic,iN), EC(ic,iN), ...
                EGF(ic,iN), EGC(ic,iN), QrelF(ic,iN), QrelC(ic,iN));
    end
end
