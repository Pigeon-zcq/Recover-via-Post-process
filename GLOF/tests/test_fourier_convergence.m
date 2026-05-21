%TEST_FOURIER_CONVERGENCE  Fourier L2/Linf convergence for x^0.5(-log x).

add_paths;

r = 0.5;
k = 1;
basis_counts = [4 8 16 32 64];
xref = linspace(0, 1, 2000);
f_handle = @(x) f_logsing(x, r, k);
fref = f_handle(xref);

err_l2 = zeros(size(basis_counts));
err_linf = zeros(size(basis_counts));

fprintf('\n[Fourier] f(x)=x^{%.1f}(-log x)^{%d}\n', r, k);
fprintf('  Nbasis   modes used                 L2 error       Linf error\n');

for iN = 1:numel(basis_counts)
    Nbasis = basis_counts(iN);
    ell = fourier_modes_for_basis_count(Nbasis);
    c = fourier_coeffs_quad(f_handle, ell);
    fN = eval_fourier_sum(xref, ell, c);

    diff = fref - fN;
    err_l2(iN) = sqrt(trapz(xref, abs(diff).^2));
    err_linf(iN) = max(abs(diff));

    fprintf('  %6d   %-24s %.6e   %.6e\n', ...
            Nbasis, sprintf('[%d,%d]', ell(1), ell(end)), ...
            err_l2(iN), err_linf(iN));
end

assert(all(isfinite(err_l2)), 'Fourier L2 errors contain non-finite values.');
assert(all(isfinite(err_linf)), 'Fourier Linf errors contain non-finite values.');
assert(err_l2(end) < err_l2(1), 'Fourier L2 error did not decrease overall.');
assert(err_linf(end) < err_linf(1), 'Fourier Linf error did not decrease overall.');

out_dir = fileparts(mfilename('fullpath'));
out_file = fullfile(out_dir, 'fourier_convergence_r05_k1.png');

fig = figure('Visible', 'off', 'Name', 'Fourier convergence r=0.5 k=1');
loglog(basis_counts, err_l2, '-o', 'LineWidth', 1.5); hold on;
loglog(basis_counts, err_linf, '-s', 'LineWidth', 1.5);
grid on;
xlabel('number of Fourier basis functions N');
ylabel('error');
legend('L2', 'Linf', 'Location', 'southwest');
title('Fourier convergence for f(x)=x^{0.5}(-log x)');
saveas(fig, out_file);
close(fig);

assert(exist(out_file, 'file') == 2, 'Fourier convergence figure was not created.');
fprintf('  Saved %s\n', out_file);

function ell = fourier_modes_for_basis_count(Nbasis)
%FOURIER_MODES_FOR_BASIS_COUNT  Return exactly Nbasis low-frequency modes.
%
%   For odd Nbasis this is the symmetric set -K:K.  For even Nbasis an
%   exactly symmetric complex-exponential set would have an odd count, so
%   the extra mode is placed on the positive side.

    left = floor((Nbasis - 1) / 2);
    right = Nbasis - left - 1;
    ell = (-left:right).';
end
