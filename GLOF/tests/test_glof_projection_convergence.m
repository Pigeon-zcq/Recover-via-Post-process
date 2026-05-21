%TEST_GLOF_PROJECTION_CONVERGENCE  Direct GLOF projection accuracy tests.
%
% This script intentionally keeps the two experiments expanded instead of
% hiding the workflow in local helper functions.  It is easier to debug by
% running one section at a time and inspecting fref, g_all, R, and the error
% arrays directly in the MATLAB workspace.

add_paths;

xref = linspace(0, 1, 2000);
basis_counts = 1:50;
alpha = 0;
lambda = 0;

%% Sanity check: stable log-domain coefficients against existing quadrature

sanity_r = 0.5;
sanity_k = 1;
sanity_m = 3;
sanity_beta = 2 * sanity_r;
sanity_Q = 30;

g_stable = glof_coeffs_logsing_stable(sanity_r, sanity_k, sanity_m, ...
                                      alpha, sanity_beta, lambda, sanity_Q);
g_quad = glof_coeffs_exact(@(x) f_logsing(x, sanity_r, sanity_k), ...
                           sanity_m, alpha, sanity_beta, lambda, sanity_Q);
sanity_rel = norm(g_stable - g_quad) / (norm(g_stable) + eps);

fprintf('\n[GLOF] stable coefficient sanity relative difference: %.3e\n', sanity_rel);
assert(sanity_rel < 1e-8, 'Stable and existing direct GLOF coefficients disagree.');

%% Experiment 1: k = 0, r varies, alpha = 0, lambda = 0, beta = 3

beta = 3;
k_fixed = 0;
r_values = [0.1 0.3 0.5 0.7 0.9];

err_l2_r = zeros(numel(r_values), numel(basis_counts));
err_linf_r = zeros(numel(r_values), numel(basis_counts));

fprintf('\n[GLOF r-sweep] k=0, alpha=0, lambda=0, beta=3\n');
fprintf('  r      bestN_L2  min L2       final L2     bestN_Linf  min Linf     final Linf\n');

for ir = 1:numel(r_values)
    r = r_values(ir);
    k = k_fixed;
    fref = f_logsing(xref, r, k);

    max_m = max(basis_counts) - 1;
    Q = max(90, max_m + k + 30);
    g_all = glof_coeffs_logsing_stable(r, k, max_m, alpha, beta, lambda, Q);

    for iN = 1:numel(basis_counts)
        Nbasis = basis_counts(iN);
        g = g_all(1:Nbasis);
        R = eval_glof_reconstruction(xref, g, alpha, beta, lambda);
        diff = fref - R;

        err_l2_r(ir, iN) = sqrt(trapz(xref, abs(diff).^2));
        err_linf_r(ir, iN) = max(abs(diff));
    end

    [best_l2, idx_l2] = min(err_l2_r(ir, :));
    [best_linf, idx_linf] = min(err_linf_r(ir, :));
    fprintf('  %.1f    %4d      %.6e  %.6e    %4d        %.6e  %.6e\n', ...
            r, basis_counts(idx_l2), best_l2, err_l2_r(ir, end), ...
            basis_counts(idx_linf), best_linf, err_linf_r(ir, end));
end

assert(all(isfinite(err_l2_r(:))), 'GLOF r-sweep L2 errors contain non-finite values.');
assert(all(isfinite(err_linf_r(:))), 'GLOF r-sweep Linf errors contain non-finite values.');

out_dir = fileparts(mfilename('fullpath'));
out_file_r = fullfile(out_dir, 'glof_projection_k0_beta3.png');

fig = figure('Visible', 'off', ...
             'Name', 'GLOF direct projection k=0 beta=3', ...
             'Position', [100 100 1000 430]);

subplot(1, 2, 1);
semilogy(basis_counts, max(err_l2_r, realmin), 'LineWidth', 1.4);
grid on;
xlabel('number of GLOF basis functions N');
ylabel('L2 error');
title('L2 convergence');
legend(compose('r=%.1f', r_values), 'Location', 'best');

subplot(1, 2, 2);
semilogy(basis_counts, max(err_linf_r, realmin), 'LineWidth', 1.4);
grid on;
xlabel('number of GLOF basis functions N');
ylabel('Linf error');
title('Linf convergence');
legend(compose('r=%.1f', r_values), 'Location', 'best');

sgtitle('GLOF direct projection, k=0, \alpha=0, \lambda=0, \beta=3');
saveas(fig, out_file_r);
close(fig);

assert(exist(out_file_r, 'file') == 2, 'GLOF r-sweep figure was not created.');
fprintf('  Saved %s\n', out_file_r);

%% Experiment 2: r = 1, k varies, alpha = 0, lambda = 0, beta = 5

beta = 5;
r_fixed = 1;
k_values = [1 2 3 4 5];

err_l2_k = zeros(numel(k_values), numel(basis_counts));
err_linf_k = zeros(numel(k_values), numel(basis_counts));

fprintf('\n[GLOF k-sweep] r=1, alpha=0, lambda=0, beta=5\n');
fprintf('  k      bestN_L2  min L2       final L2     bestN_Linf  min Linf     final Linf\n');

for ik = 1:numel(k_values)
    r = r_fixed;
    k = k_values(ik);
    fref = f_logsing(xref, r, k);

    max_m = max(basis_counts) - 1;
    Q = max(90, max_m + k + 30);
    g_all = glof_coeffs_logsing_stable(r, k, max_m, alpha, beta, lambda, Q);

    for iN = 1:numel(basis_counts)
        Nbasis = basis_counts(iN);
        g = g_all(1:Nbasis);
        R = eval_glof_reconstruction(xref, g, alpha, beta, lambda);
        diff = fref - R;

        err_l2_k(ik, iN) = sqrt(trapz(xref, abs(diff).^2));
        err_linf_k(ik, iN) = max(abs(diff));
    end

    [best_l2, idx_l2] = min(err_l2_k(ik, :));
    [best_linf, idx_linf] = min(err_linf_k(ik, :));
    fprintf('  %d      %4d      %.6e  %.6e    %4d        %.6e  %.6e\n', ...
            k, basis_counts(idx_l2), best_l2, err_l2_k(ik, end), ...
            basis_counts(idx_linf), best_linf, err_linf_k(ik, end));
end

assert(all(isfinite(err_l2_k(:))), 'GLOF k-sweep L2 errors contain non-finite values.');
assert(all(isfinite(err_linf_k(:))), 'GLOF k-sweep Linf errors contain non-finite values.');

out_file_k = fullfile(out_dir, 'glof_projection_r1_beta5.png');

fig = figure('Visible', 'off', ...
             'Name', 'GLOF direct projection r=1 beta=5', ...
             'Position', [100 100 1000 430]);

subplot(1, 2, 1);
semilogy(basis_counts, max(err_l2_k, realmin), 'LineWidth', 1.4);
grid on;
xlabel('number of GLOF basis functions N');
ylabel('L2 error');
title('L2 convergence');
legend(compose('k=%d', k_values), 'Location', 'best');

subplot(1, 2, 2);
semilogy(basis_counts, max(err_linf_k, realmin), 'LineWidth', 1.4);
grid on;
xlabel('number of GLOF basis functions N');
ylabel('Linf error');
title('Linf convergence');
legend(compose('k=%d', k_values), 'Location', 'best');

sgtitle('GLOF direct projection, r=1, \alpha=0, \lambda=0, \beta=5');
saveas(fig, out_file_k);
close(fig);

assert(exist(out_file_k, 'file') == 2, 'GLOF k-sweep figure was not created.');
fprintf('  Saved %s\n', out_file_k);
