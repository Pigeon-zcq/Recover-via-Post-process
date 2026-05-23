% diag4_ex65.m — test non-integer lambda and larger N for Example 6.5
addpath('/Users/zcq/PaperCode/Recover/Gegenbauer');
clear; clc;

p=1; q=2; bet=0.5; a_sub=-0.5; b_sub=0.0; Tfinal=log(2);
f_et=@(x) sqrt(max(0,1+2*x));

fprintf('=== Test A: Exact non-integer lambda = N/80 for large N ===\n');
fprintf('  %4s  %8s  %4s  %10s\n', 'N', 'lambda', 'm', 'maxerr');
for N = [40, 80, 160, 320]
    lam = N/80;  % non-integer for N<80
    m   = floor(N/4);
    nt  = 3*N; nzn = N;
    hf = FourierGalerkinTransport(N, p, q, Tfinal, 0.5);
    gt = compute_gt(lam, max(m,nt)+1);
    hg = gcoe(hf, N, m, lam, bet, a_sub, b_sub, gt, nt);
    [u, xx] = gegenbauer_reconstruct(hg, m, lam, bet, a_sub, b_sub, nzn, gt);
    ex = arrayfun(f_et, xx);
    fprintf('  %4d  %8.4f  %4d  %10.3e\n', N, lam, m, max(abs(u-ex)));
end

fprintf('\n=== Test B: Fixed lambda=1 but larger N (check convergence) ===\n');
fprintf('  %4s  %8s  %4s  %10s\n', 'N', 'lambda', 'm', 'maxerr');
for N = [40, 80, 160, 320]
    lam = 1;
    m   = floor(N/4);
    nt  = 3*N; nzn = N;
    hf = FourierGalerkinTransport(N, p, q, Tfinal, 0.5);
    gt = compute_gt(lam, max(m,nt)+1);
    hg = gcoe(hf, N, m, lam, bet, a_sub, b_sub, gt, nt);
    [u, xx] = gegenbauer_reconstruct(hg, m, lam, bet, a_sub, b_sub, nzn, gt);
    ex = arrayfun(f_et, xx);
    fprintf('  %4d  %8.4f  %4d  %10.3e\n', N, lam, m, max(abs(u-ex)));
end

fprintf('\n=== Test C: Exact hf, non-integer lambda = N/80 ===\n');
u_T = @(x) piecewise_u_T(x);
fprintf('  %4s  %8s  %4s  %10s\n', 'N', 'lambda', 'm', 'maxerr');
for N = [40, 80, 160]
    lam = N/80;
    m   = floor(N/4);
    nt  = 3*N; nzn = N;
    hf = getPointValuefromFC(N, u_T);
    gt = compute_gt(lam, max(m,nt)+1);
    hg = gcoe(hf, N, m, lam, bet, a_sub, b_sub, gt, nt);
    [u, xx] = gegenbauer_reconstruct(hg, m, lam, bet, a_sub, b_sub, nzn, gt);
    ex = arrayfun(f_et, xx);
    fprintf('  %4d  %8.4f  %4d  %10.3e\n', N, lam, m, max(abs(u-ex)));
end

function v = piecewise_u_T(x)
    if x >= 0.5
        v = sqrt(2*x-1);
    elseif x >= -0.5
        v = sqrt(max(0,1+2*x));
    else
        v = sqrt(max(0,3+2*x));
    end
end
