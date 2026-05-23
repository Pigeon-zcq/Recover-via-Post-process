% test_GegenbaurReconstruct.m
%
% Main test script for the Gegenbauer reconstruction algorithm.
% Reproduces Tables 1, 2 and Example 6.5 from:
%   Chen & Shu (2015), J. Sci. Comput. 65, 1145-1165
%
% Required files in the same directory:
%   compute_gt.m          getPointValuefromFC.m
%   transf_y2x.m          transf_x2y.m
%   gegen.m               gcoe.m
%   gegenbauer_reconstruct.m
%   run_reconstruction.m  findBestParameters.m
%   FourierGalerkinTransport.m
%   f_example1.m          f_transport_ic.m      f_transport_exact.m
%
% Run with:  >> test_GegenbaurReconstruct

clear; clc;
fprintf('=== Gegenbauer Reconstruction (Chen & Shu 2015) ===\n\n');

% -----------------------------------------------------------------------
%  PART 1 – Example 6.1: u(x) = cos(x) + sin(x)*sqrt(1+x),  p=1, q=2
%  Linear parameter choice: lambda = N/16, m = 3N/80  (Table 1)
% -----------------------------------------------------------------------
fprintf('--- Example 6.1: u(x) = cos(x) + sin(x)*sqrt(1+x) ---\n');
fprintf('  Linear parameter choice (lambda = N/16, m = 3N/80)\n');
p = 1;  q = 2;
f1 = @(x) f_example1(x, p, q);

fprintf('  %6s  %6s  %4s  %10s\n', 'N', 'lambda', 'm', 'L_inf err');
Nlist = [40, 80, 160, 320];
for N = Nlist
    lam = max(1, floor(N/16));
    m   = max(1, floor(3*N/80));
    err = run_reconstruction(N, lam, m, p, q, -1, 1, f1, f1);
    fprintf('  %6d  %6d  %4d  %10.3e\n', N, lam, m, err);
end

% -----------------------------------------------------------------------
%  PART 2 – Example 6.1: optimal (lambda, m) search  (N = 40, 80)
% -----------------------------------------------------------------------
fprintf('\n  Optimal parameter search:\n');
fprintf('  %6s  %6s  %4s  %10s\n', 'N', 'lambda', 'm', 'L_inf err');
for N = [40, 80]
    [mb, lb, eb] = findBestParameters(N, p, q, -1, 1, f1, f1);
    fprintf('  %6d  %6d  %4d  %10.3e\n', N, lb, mb, eb);
end

% -----------------------------------------------------------------------
%  PART 3 – Example 6.2: u(x) = cos(x) + sin(x)*(1+x)^(1/3),  p=1, q=3
%  Linear parameter choice: lambda = N/8, m = N/32  (Table 2)
% -----------------------------------------------------------------------
fprintf('\n--- Example 6.2: u(x) = cos(x) + sin(x)*(1+x)^(1/3) ---\n');
fprintf('  Linear parameter choice (lambda = N/8, m = N/32)\n');
p = 1;  q = 3;
f2 = @(x) f_example1(x, p, q);

fprintf('  %6s  %6s  %4s  %10s\n', 'N', 'lambda', 'm', 'L_inf err');
for N = Nlist
    lam = max(1, floor(N/8));
    m   = max(1, floor(N/32));
    err = run_reconstruction(N, lam, m, p, q, -1, 1, f2, f2);
    fprintf('  %6d  %6d  %4d  %10.3e\n', N, lam, m, err);
end

% -----------------------------------------------------------------------
%  PART 4 – Example 6.5: transport PDE  u_t - x*u_x = 0
%  Recover u(x, log2) = sqrt(1+2x) on [-0.5, 0]  (Table 7)
% -----------------------------------------------------------------------
fprintf('\n--- Example 6.5: Transport PDE u_t - xu_x = 0 ---\n');
fprintf('  Recovering u(x,log2)=sqrt(1+2x) on [-0.5,0]\n');
p = 1;  q = 2;
Tfinal = log(2);
a_sub  = -0.5;
b_sub  =  0.0;
bet    = 1/q;
f_et   = @(x) f_transport_exact(x, p, q);

fprintf('  %6s  %6s  %4s  %10s\n', 'N', 'lambda', 'm', 'L_inf err');
for N = [10, 20, 40, 80]
    % Linear choice from Table 7: lambda = N/80, m = N/4
    lam = max(1, floor(N/80));
    m   = max(1, floor(N/4));
    nt  = 3*N;
    nzn = N;

    hf      = FourierGalerkinTransport(N, p, q, Tfinal, 0.5);
    gt      = compute_gt(lam, max(m, nt) + 1);
    hg      = gcoe(hf, N, m, lam, bet, a_sub, b_sub, gt, nt);
    [u, xx] = gegenbauer_reconstruct(hg, m, lam, bet, a_sub, b_sub, nzn, gt);
    ex      = arrayfun(f_et, xx);
    fprintf('  %6d  %6d  %4d  %10.3e\n', N, lam, m, max(abs(u - ex)));
end

% -----------------------------------------------------------------------
%  PART 5 – Pointwise error plot for Example 6.1 with N=80
% -----------------------------------------------------------------------
fprintf('\nGenerating pointwise error plot (Example 6.1, N=80)...\n');
p = 1;  q = 2;
f1 = @(x) f_example1(x, p, q);
N  = 80;
lam = max(1, floor(N/16));
m   = max(1, floor(3*N/80));
bet = 1/q;
nt  = 3*N;
nzn = N;

hf      = getPointValuefromFC(N, f1);
gt      = compute_gt(lam, max(m, nt) + 1);
hg      = gcoe(hf, N, m, lam, bet, -1, 1, gt, nt);
[u, xx] = gegenbauer_reconstruct(hg, m, lam, bet, -1, 1, nzn, gt);
ex      = arrayfun(f1, xx);

figure;
semilogy(xx, abs(u - ex) + eps, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2);
xlabel('x');
ylabel('|error|');
title(sprintf('Pointwise error: Example 6.1  (N=%d, \\lambda=%d, m=%d)', N, lam, m));
grid on;

fprintf('\nDone.\n');
