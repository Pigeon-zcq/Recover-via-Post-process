%TEST_GLOF_FROM_FOURIER_SIN_MULTIPLIER  GLOF-F test for sin(x)x^r(-log x)^k.

fprintf('\n[GLOF from Fourier: sin multiplier] f(x)=sin(x)x^r(-log x)^k\n');

r = 0.5;
k = 1;
gamma_par = 0.10;

Ns = [16 32 64 128 256];
xref = linspace(0, 1, 2000);
f_handle = @(x) sin(x) .* f_logsing(x, r, k);
fref = f_handle(xref);

theta_values = [r, r + 1];
theta_labels = {'theta = r', 'theta = r+1'};

err_fourier_l2 = zeros(size(Ns));
err_fourier_linf = zeros(size(Ns));
err_glof_l2 = zeros(numel(theta_values), numel(Ns));
err_glof_linf = zeros(numel(theta_values), numel(Ns));
glof_basis_counts = zeros(size(Ns));
Qs = zeros(size(Ns));

for iN = 1:numel(Ns)
    N = Ns(iN);
    ell = (-N:N).';
    c = fourier_coeffs_quad(f_handle, ell, 1e-11, 1e-13);

    f_fourier = eval_fourier_sum(xref, ell, c);
    diff_fourier = f_fourier - fref;
    err_fourier_l2(iN) = sqrt(trapz(xref, abs(diff_fourier).^2));
    err_fourier_linf(iN) = max(abs(diff_fourier));

    M = max(k + 1, ceil(k + gamma_par * N));
    m = M - 1;
    Q = max(M + 4, ceil(0.10 * N));
    glof_basis_counts(iN) = M;
    Qs(iN) = Q;

    for itheta = 1:numel(theta_values)
        theta = theta_values(itheta);
        alpha = 0;
        lambda = 0;
        beta = 2 * theta;

        g = glof_coeffs_from_fourier(ell, c, m, alpha, beta, lambda, Q);
        R = eval_glof_reconstruction(xref, g, alpha, beta, lambda);
        diff_glof = R - fref;

        err_glof_l2(itheta, iN) = sqrt(trapz(xref, abs(diff_glof).^2));
        err_glof_linf(itheta, iN) = max(abs(diff_glof));
    end

    fprintf(['  N=%3d modes=%4d M=%2d Q=%2d | Fourier L2 %.3e Linf %.3e | ' ...
             'GLOF theta=r Linf %.3e | GLOF theta=r+1 Linf %.3e\n'], ...
            N, 2*N+1, M, Q, err_fourier_l2(iN), err_fourier_linf(iN), ...
            err_glof_linf(1, iN), err_glof_linf(2, iN));
end

assert(all(isfinite(err_fourier_l2)), 'Fourier L2 errors contain non-finite values.');
assert(all(isfinite(err_fourier_linf)), 'Fourier Linf errors contain non-finite values.');
assert(all(isfinite(err_glof_l2(:))), 'GLOF-F L2 errors contain non-finite values.');
assert(all(isfinite(err_glof_linf(:))), 'GLOF-F Linf errors contain non-finite values.');

assert(err_glof_linf(1, end) < err_glof_linf(1, 1), ...
       'GLOF-F theta=r Linf error did not decrease overall.');
assert(err_glof_linf(2, end) < err_glof_linf(2, 1), ...
       'GLOF-F theta=r+1 Linf error did not decrease overall.');
assert(err_glof_l2(1, end) < err_glof_l2(1, 1), ...
       'GLOF-F theta=r L2 error did not decrease overall.');
assert(err_glof_l2(2, end) < err_glof_l2(2, 1), ...
       'GLOF-F theta=r+1 L2 error did not decrease overall.');

out_dir = fileparts(mfilename('fullpath'));
out_file = fullfile(out_dir, 'glof_from_fourier_sin_multiplier_r05_k1.png');

fig = figure('Visible', 'off', ...
             'Name', 'GLOF-F sin multiplier r=0.5 k=1', ...
             'Position', [100 100 1100 450]);

subplot(1, 2, 1);
semilogy(Ns, err_fourier_l2, '-o', 'LineWidth', 1.4); hold on;
semilogy(Ns, err_glof_l2(1, :), '-^', 'LineWidth', 1.4);
semilogy(Ns, err_glof_l2(2, :), '-d', 'LineWidth', 1.4);
grid on;
xlabel('Fourier half-bandwidth N');
ylabel('L2 error on [0,1]');
title('L2 reconstruction error');
legend('Fourier', theta_labels{1}, theta_labels{2}, 'Location', 'southwest');

subplot(1, 2, 2);
semilogy(Ns, err_fourier_linf, '-o', 'LineWidth', 1.4); hold on;
semilogy(Ns, err_glof_linf(1, :), '-^', 'LineWidth', 1.4);
semilogy(Ns, err_glof_linf(2, :), '-d', 'LineWidth', 1.4);
grid on;
xlabel('Fourier half-bandwidth N');
ylabel('Linf error on [0,1]');
title('Linf reconstruction error');
legend('Fourier', theta_labels{1}, theta_labels{2}, 'Location', 'southwest');

sgtitle(sprintf('GLOF-F for sin(x)x^{%.1f}(-log x)^{%d}, M=ceil(k+%.2fN)', ...
                r, k, gamma_par));
print(fig, out_file, '-dpng', '-r200');
close(fig);

assert(exist(out_file, 'file') == 2, ...
       'GLOF-F sin multiplier figure was not created.');

fprintf('[GLOF from Fourier: sin multiplier] saved %s\n', out_file);
fprintf('[GLOF from Fourier: sin multiplier] GLOF basis counts: %s\n', ...
        mat2str(glof_basis_counts));
