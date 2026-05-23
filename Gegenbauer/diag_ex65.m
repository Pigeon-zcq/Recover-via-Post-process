% diag_ex65.m — diagnose Example 6.5
% Goal: isolate whether the bug is in FourierGalerkinTransport or in gcoe/reconstruct

addpath('/Users/zcq/PaperCode/Recover/Gegenbauer');
clear; clc;

p = 1; q = 2; bet = 0.5;
a_sub = -0.5; b_sub = 0.0;
Tfinal = log(2);

% --- Exact solution at T=log(2) on the FULL periodic domain [-1,1] ---
% u_T(x) = (3+2x)^.5  x in [-1,-0.5)
%          (1+2x)^.5  x in [-0.5,0.5)
%          (2x-1)^.5  x in [0.5,1)
u_T = @(x) exact_sol_T(x);

f_et = @(x) max(0, 1+2*x)^0.5;

fprintf('=== DIAGNOSTIC: Example 6.5 ===\n\n');
fprintf('--- Test A: Use exact Fourier coefficients of u_T(x) ---\n');
fprintf('  (If this works, bug is in FourierGalerkinTransport)\n\n');
fprintf('  %4s  %6s  %4s  %10s  %10s\n', 'N', 'lambda', 'm', 'err(exact hf)', 'err(galerkin)');

for N = [10, 20, 40]
    lam = max(1, floor(N/80));
    m   = max(1, floor(N/4));
    nt  = 3*N;
    nzn = N;

    % Exact Fourier coefficients of u_T
    hf_exact = getPointValuefromFC(N, u_T);

    % Galerkin Fourier coefficients
    hf_gal = FourierGalerkinTransport(N, p, q, Tfinal, 0.5);

    gt = compute_gt(lam, max(m, nt) + 1);

    % Reconstruction from exact hf
    hg_e  = gcoe(hf_exact, N, m, lam, bet, a_sub, b_sub, gt, nt);
    [u_e, xx_e] = gegenbauer_reconstruct(hg_e, m, lam, bet, a_sub, b_sub, nzn, gt);
    ex_e  = arrayfun(f_et, xx_e);
    err_e = max(abs(u_e - ex_e));

    % Reconstruction from Galerkin hf
    hg_g  = gcoe(hf_gal, N, m, lam, bet, a_sub, b_sub, gt, nt);
    [u_g, xx_g] = gegenbauer_reconstruct(hg_g, m, lam, bet, a_sub, b_sub, nzn, gt);
    ex_g  = arrayfun(f_et, xx_g);
    err_g = max(abs(u_g - ex_g));

    fprintf('  %4d  %6d  %4d  %10.3e  %10.3e\n', N, lam, m, err_e, err_g);
end

fprintf('\n--- Test B: Fourier coefficients of k=0..N of hf_galerkin ---\n');
N = 10;
hf_gal = FourierGalerkinTransport(N, p, q, Tfinal, 0.5);
hf_exact = getPointValuefromFC(N, u_T);
hf_ic    = getPointValuefromFC(N, @(x) max(0,1+x)^0.5);

fprintf('  k  |hf_gal|    |hf_exact|   |hf_ic|\n');
for k = 0:min(N,5)
    fprintf('  %2d  %.4e   %.4e   %.4e\n', k, ...
        abs(hf_gal(k+N+1)), abs(hf_exact(k+N+1)), abs(hf_ic(k+N+1)));
end

fprintf('\n--- Test C: Direct evaluation of Fourier series at test points ---\n');
x_test = [-0.4, -0.3, -0.2, -0.1, 0.0];
fprintf('  x       exact      f_N(x,gal)   f_N(x,exact_hf)\n');
for x = x_test
    fN_gal   = 0; fN_ex = 0;
    for k = -N:N
        fN_gal = fN_gal + real(hf_gal(k+N+1)   * exp(1i*k*pi*x));
        fN_ex  = fN_ex  + real(hf_exact(k+N+1)  * exp(1i*k*pi*x));
    end
    fprintf('  %5.2f  %9.5f  %11.5f  %11.5f\n', x, f_et(x), fN_gal, fN_ex);
end

function v = exact_sol_T(x)
    if x >= 0.5
        v = sqrt(2*x-1);
    elseif x >= -0.5
        v = sqrt(1+2*x);
    else
        v = sqrt(3+2*x);
    end
end
