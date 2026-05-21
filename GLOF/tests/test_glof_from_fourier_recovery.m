%TEST_GLOF_FROM_FOURIER_RECOVERY  Recover GLOF accuracy from Fourier data.
%
% This test exercises the core post-processing pipeline
%
%   Fourier coefficients -> GLOF coefficients -> reconstructed function.
%
% The reconstruction coefficients are computed only by
% GLOF_COEFFS_FROM_FOURIER.  The stable direct projection is used below only
% as a coefficient diagnostic for the known model function.

fprintf('\n[GLOF from Fourier] recovering f(x)=x^r(-log x)^k\n');

r = 0.5;
k = 1;

% N is the Fourier half-bandwidth, so each run uses 2*N+1 Fourier modes.
Ns = [16 32 64 128 256];

% Chen-Shen style N-dependent parameters.
a_par = 0.08;
b_par = 0.20;
q_factor = 1.50;
gamma_par = 0.10;

xref = linspace(0, 1, 2000);
f_handle = @(x) f_logsing(x, r, k);
fref = f_handle(xref);

err_fourier_l2 = zeros(size(Ns));
err_fourier_linf = zeros(size(Ns));
err_glof_l2 = zeros(size(Ns));
err_glof_linf = zeros(size(Ns));
coeff_rel = zeros(size(Ns));
q_rel = zeros(size(Ns));
glof_basis_counts = zeros(size(Ns));
Qs = zeros(size(Ns));

for iN = 1:numel(Ns)
    N = Ns(iN);

    alpha = a_par * N;
    beta = 2*r + b_par * N;
    lambda = b_par * N;

    % In the core routines, m is the highest GLOF degree and the number of
    % basis functions is m+1.  The experiment prescription uses the basis
    % count M=ceil(k+gamma*N), so the core degree is M-1.
    M = ceil(k + gamma_par * N);
    m = M - 1;
    Q = ceil(q_factor * N);

    glof_basis_counts(iN) = M;
    Qs(iN) = Q;

    ell = (-N:N).';
    c = fourier_coeffs_quad(f_handle, ell);

    f_fourier = eval_fourier_sum(xref, ell, c);
    diff_fourier = f_fourier - fref;
    err_fourier_l2(iN) = sqrt(trapz(xref, abs(diff_fourier).^2));
    err_fourier_linf(iN) = max(abs(diff_fourier));

    g_fourier = glof_coeffs_from_fourier(ell, c, m, alpha, beta, lambda, Q);
    g_fourier_2Q = glof_coeffs_from_fourier(ell, c, m, alpha, beta, lambda, 2*Q);
    g_direct = glof_coeffs_logsing_stable(r, k, m, alpha, beta, lambda, 2*Q);

    R = eval_glof_reconstruction(xref, g_fourier, alpha, beta, lambda);
    diff_glof = R - fref;
    err_glof_l2(iN) = sqrt(trapz(xref, abs(diff_glof).^2));
    err_glof_linf(iN) = max(abs(diff_glof));

    coeff_rel(iN) = norm(g_fourier - g_direct) / (norm(g_direct) + 1e-300);
    q_rel(iN) = norm(g_fourier - g_fourier_2Q) / (norm(g_fourier_2Q) + 1e-300);

    fprintf(['  N=%3d modes=%4d | alpha=%.2f beta=%.2f lambda=%.2f M=%d Q=%d | ' ...
             'Fourier L2 %.3e Linf %.3e | GLOF-F L2 %.3e Linf %.3e | ' ...
             'coeff %.3e Qrel %.3e\n'], ...
            N, 2*N+1, alpha, beta, lambda, M, Q, ...
            err_fourier_l2(iN), err_fourier_linf(iN), ...
            err_glof_l2(iN), err_glof_linf(iN), coeff_rel(iN), q_rel(iN));
end

assert(all(isfinite(err_fourier_l2)), 'Fourier L2 errors contain non-finite values.');
assert(all(isfinite(err_fourier_linf)), 'Fourier Linf errors contain non-finite values.');
assert(all(isfinite(err_glof_l2)), 'GLOF-from-Fourier L2 errors contain non-finite values.');
assert(all(isfinite(err_glof_linf)), 'GLOF-from-Fourier Linf errors contain non-finite values.');
assert(all(isfinite(coeff_rel)), 'Coefficient diagnostics contain non-finite values.');
assert(all(isfinite(q_rel)), 'Q-convergence diagnostics contain non-finite values.');

% Diagnose whether the post-processed result recovers far more accuracy than
% the raw Fourier partial sum and whether the high-N tail has exponential
% type decay.  With M=ceil(k+0.1*N), this test intentionally exposes the
% high-mode endpoint-noise amplification if it appears.
tail_fit_linf = polyfit(Ns(end-2:end), log(max(err_glof_linf(end-2:end), realmin)), 1);
tail_fit_l2 = polyfit(Ns(end-2:end), log(max(err_glof_l2(end-2:end), realmin)), 1);

recovers_linf = err_glof_linf(end) < 1e-3 * err_fourier_linf(end);
recovers_l2 = err_glof_l2(end) < 1e-3 * err_fourier_l2(end);
tail_exponential = (tail_fit_linf(1) < -2e-2) && (tail_fit_l2(1) < -2e-2);

if ~(recovers_linf && recovers_l2 && tail_exponential)
    warning('GLOF:FourierRecoveryNotExponential', ...
            ['M=ceil(k+%.2f*N) did not recover exponential accuracy for this ' ...
             'core Fourier-to-GLOF pipeline. Final ratios: Linf %.3e, L2 %.3e. ' ...
             'Tail slopes: Linf %.3e, L2 %.3e.'], ...
            gamma_par, ...
            err_glof_linf(end)/(err_fourier_linf(end) + 1e-300), ...
            err_glof_l2(end)/(err_fourier_l2(end) + 1e-300), ...
            tail_fit_linf(1), tail_fit_l2(1));
end

out_dir = fileparts(mfilename('fullpath'));
out_file = fullfile(out_dir, 'glof_from_fourier_recovery_r05_k1.png');

fig = figure('Visible', 'off', ...
             'Name', 'GLOF recovery from Fourier data r=0.5 k=1', ...
             'Position', [100 100 1100 450]);

subplot(1, 2, 1);
semilogy(Ns, err_fourier_l2, '-o', 'LineWidth', 1.4); hold on;
semilogy(Ns, err_fourier_linf, '-s', 'LineWidth', 1.4);
semilogy(Ns, err_glof_l2, '-^', 'LineWidth', 1.4);
semilogy(Ns, err_glof_linf, '-d', 'LineWidth', 1.4);
grid on;
xlabel('Fourier half-bandwidth N');
ylabel('error on [0,1]');
title('Function reconstruction error');
legend('Fourier L2', 'Fourier Linf', 'GLOF-F L2', 'GLOF-F Linf', ...
       'Location', 'southwest');

subplot(1, 2, 2);
semilogy(Ns, coeff_rel, '-o', 'LineWidth', 1.4); hold on;
semilogy(Ns, q_rel, '-s', 'LineWidth', 1.4);
grid on;
xlabel('Fourier half-bandwidth N');
ylabel('relative coefficient error');
title('GLOF coefficient diagnostics');
legend('vs direct projection', 'Q vs 2Q', 'Location', 'southwest');

sgtitle(sprintf(['Fourier coefficients to GLOF recovery, r=%.1f, k=%d, ' ...
                 'a=%.2f, b=%.2f, M=ceil(k+%.2fN)'], ...
                r, k, a_par, b_par, gamma_par));
print(fig, out_file, '-dpng', '-r200');
close(fig);

assert(exist(out_file, 'file') == 2, ...
       'GLOF-from-Fourier recovery figure was not created.');

fprintf('[GLOF from Fourier] saved %s\n', out_file);
fprintf('[GLOF from Fourier] tail slopes: Linf %.3e, L2 %.3e\n', ...
        tail_fit_linf(1), tail_fit_l2(1));
fprintf('[GLOF from Fourier] GLOF basis counts: %s\n', ...
        mat2str(glof_basis_counts));
