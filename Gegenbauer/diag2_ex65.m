% diag2_ex65.m — minimal diagnostic for gcoe / gegenbauer_reconstruct
% Tests the reconstruction step independently of FourierGalerkinTransport

addpath('/Users/zcq/PaperCode/Recover/Gegenbauer');
clear; clc;

bet = 0.5; a = -0.5; b = 0.0;
lam = 1; m = 4; nt = 30; nzn = 10;
N = 10;

% --- exact function F(y) = (1+y)/2 after coordinate transform
F_exact = @(y) (1+y)/2;

% --- Analytic Gegenbauer coefficients of F(y) = (1+y)/2 with lambda=1:
% hg(0) = 0.5, hg(1) = 0.25, hg(l>=2) = 0
hg_analytic = [0.5, 0.25, zeros(1, m-1)];
fprintf('Analytic hg: '); fprintf('%.4f ', hg_analytic); fprintf('\n');

% --- Test 1: gegenbauer_reconstruct with analytic hg (bypass gcoe)
gt = compute_gt(lam, max(m, nt) + 1);
[u_ana, xx_ana] = gegenbauer_reconstruct(hg_analytic, m, lam, bet, a, b, nzn, gt);
ex_ana = arrayfun(@(x) sqrt(max(0,1+2*x)), xx_ana);
fprintf('\nTest 1: gegenbauer_reconstruct with ANALYTIC hg\n');
fprintf('  max error = %.4e\n', max(abs(u_ana - ex_ana)));
fprintf('  (expected ~machine epsilon)\n\n');

% --- Test 2: gcoe with EXACT function values (not Fourier series)
xi = cos(pi * (1:nt-1) / nt);
xr = transf_y2x(xi, bet, a, b);
ff_exact = arrayfun(F_exact, xi);  % exact F(y) at quadrature nodes

% manually apply gcoe formula
log_hlm = log(pi) - 2*gammaln(lam) - (2*lam-1)*log(2);
hg_from_exact = zeros(1, m+1);
for l = 0:m
    [~, cnl_all] = gegen(l, nt, gt);
    cnl = cnl_all(2:nt);
    fi = cnl .* ff_exact .* (1-xi).^lam .* (1+xi).^lam;
    log_hll = log_hlm + gammaln(l+2*lam) - log(l+lam) - gammaln(l+1);
    hg_from_exact(l+1) = sum(fi) * (pi/nt) * exp(-log_hll);
end
fprintf('Test 2: gcoe with EXACT F(y) values at quadrature nodes\n');
fprintf('  hg = '); fprintf('%.4f ', hg_from_exact); fprintf('\n');
fprintf('  (expected [0.5 0.25 0 0 0])\n\n');

% --- Test 3: gcoe with FOURIER-SERIES-evaluated values
hf_exact = getPointValuefromFC(N, @(x) u_T_piecewise(x));
ff_fourier = zeros(1, nt-1);
for k = -N:N
    ff_fourier = ff_fourier + real(hf_exact(k+N+1) * exp(1i*k*pi*xr));
end
hg_from_fourier = zeros(1, m+1);
for l = 0:m
    [~, cnl_all] = gegen(l, nt, gt);
    cnl = cnl_all(2:nt);
    fi = cnl .* ff_fourier .* (1-xi).^lam .* (1+xi).^lam;
    log_hll = log_hlm + gammaln(l+2*lam) - log(l+lam) - gammaln(l+1);
    hg_from_fourier(l+1) = sum(fi) * (pi/nt) * exp(-log_hll);
end
fprintf('Test 3: gcoe with FOURIER SERIES values at quadrature nodes in x-space\n');
fprintf('  hg = '); fprintf('%.4f ', hg_from_fourier); fprintf('\n');
fprintf('  (expected near [0.5 0.25 0 0 0])\n\n');

% --- Show actual ff values (Fourier series vs exact F(y)) at quad nodes ---
fprintf('Test 4: f_N(xr) vs F(y) at first 5 quadrature nodes\n');
fprintf('  i      y=xi       xr       F(y)=exact  f_N(xr)\n');
for i = 1:5
    fprintf('  %2d  %7.4f  %7.4f  %9.4f  %9.4f\n', ...
        i, xi(i), xr(i), ff_exact(i), ff_fourier(i));
end

% --- Test 5: Full gcoe call ---
fprintf('\nTest 5: Full gcoe call\n');
hg_gcoe = gcoe(hf_exact, N, m, lam, bet, a, b, gt, nt);
fprintf('  hg from gcoe = '); fprintf('%.4f ', hg_gcoe); fprintf('\n');
[u_gcoe, xx_gcoe] = gegenbauer_reconstruct(hg_gcoe, m, lam, bet, a, b, nzn, gt);
ex_gcoe = arrayfun(@(x) sqrt(max(0,1+2*x)), xx_gcoe);
fprintf('  max error from gcoe = %.4e\n', max(abs(u_gcoe - ex_gcoe)));

function v = u_T_piecewise(x)
    if x >= 0.5
        v = sqrt(2*x-1);
    elseif x >= -0.5
        v = sqrt(max(0, 1+2*x));
    else
        v = sqrt(max(0, 3+2*x));
    end
end
