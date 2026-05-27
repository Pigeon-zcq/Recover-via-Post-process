%% main_glof_reconstruct.m
% Fourier-GLOF reconstruction — script form for debugging.
%
% Reconstructs  f(x) = x^r * (-log x)^k  from finite Fourier coefficients
% using the GLOF (Generalized Laguerre-based Orthogonal Functions) method.
%
% Reference: algorithm document code.md, Sections 1–17.

clear; clc; close all;

%% ================================================================
%  Section 0 — User parameters
% ================================================================
r     = 0.5;   % power singularity exponent (r > 0)
k     = 1;     % log power  (k = 0: pure power, k >= 1: log factor)
alpha = 0.5;   % Laguerre weight parameter  (alpha > -1)
NF    = 32;    % Fourier half-bandwidth

theta       = 0.05;   % lambda = floor(theta * NF)
gamma_coeff = 0.05;   % M      = floor(gamma_coeff * NF)
Q_factor    = 4;      % Q      = Q_factor * NF  (quadrature oversampling)

%% ================================================================
%  Derived parameters  (Section 2 / Section 14.1)
% ================================================================
lambda = floor(theta * NF);    % GLOF weight exponent
beta   = lambda + 2*r;         % log-mapping param; matches x^r singularity
M      = floor(gamma_coeff * NF);
Q      = Q_factor * NF;
nu     = r;                    % = (beta - lambda)/2 under singularity matching

%% --- Parameter sanity check (Step 1) ---
assert(beta > -1,         'Need beta > -1  (got %.4f)', beta);
assert(2*r + lambda > -1, 'Need 2r+lambda > -1  (got %.4f)', 2*r+lambda);
assert(M <= Q,            'Need M <= Q  (%d vs %d)', M, Q);

fprintf('=== Fourier-GLOF Reconstruction ===\n');
fprintf('f(x) = x^%.2f * (-log x)^%d\n', r, k);
fprintf('NF=%d,  M=%d,  Q=%d\n', NF, M, Q);
fprintf('alpha=%.2f,  lambda=%d,  beta=%.4f,  nu=%.4f\n\n', alpha, lambda, beta, nu);

%% ================================================================
%  Exact function & evaluation grid
% ================================================================
f_exact = @(x) x.^r .* (-log(x)).^k;

z_eval = (0.005 : 0.005 : 0.995)';   % 199 test points, avoid x=0 and x=1
f_true = f_exact(z_eval);

%% ================================================================
%  Step 0 — Compute Fourier coefficients by adaptive quadrature
% ================================================================
fprintf('[Step 0] Computing %d Fourier coefficients ...\n', 2*NF+1);
n_modes = (-NF : NF)';
f_hat   = zeros(2*NF+1, 1);

for idx = 1 : length(n_modes)
    n = n_modes(idx);
    re_int = @(x) x.^r .* (-log(x)).^k .* cos(2*pi*n*x);
    im_int = @(x) x.^r .* (-log(x)).^k .* sin(2*pi*n*x);
    % Use Waypoints near 0 to help adaptive integrator with singularity
    wpts = [1e-10, 1e-6, 1e-3];
    I_re = integral(re_int, 0, 1, 'RelTol',1e-12, 'AbsTol',1e-14, 'Waypoints',wpts);
    I_im = integral(im_int, 0, 1, 'RelTol',1e-12, 'AbsTol',1e-14, 'Waypoints',wpts);
    f_hat(idx) = I_re - 1i*I_im;
end
fprintf('  f_hat(0) = %.8f  (exact: %.8f)\n', real(f_hat(NF+1)), ...
    integral(@(x) x.^r .* (-log(x)).^k, 0, 1));

%% ================================================================
%  Step 2 — Laguerre-Gauss nodes and weights
% ================================================================
fprintf('[Step 2] Laguerre-Gauss nodes (%d pts) ...\n', Q+1);
[y_j, omega_j] = laguerre_gauss_nodes(Q, alpha);

%% ================================================================
%  Step 3 — Map Laguerre nodes to x in (0,1)
% ================================================================
x_j = exp(-y_j / (beta + 1));

fprintf('  x_j range: [%.2e, %.4f]\n', min(x_j), max(x_j));

%% ================================================================
%  Step 4 — GLOF quadrature weights
% ================================================================
chi_lof = (beta + 1)^(-(alpha+1)) * omega_j;     % LOF weights
chi_j   = x_j.^(lambda - beta) .* chi_lof;       % GLOF weights (Section 4.4)

%% ================================================================
%  Step 5 — Fourier partial sum at quadrature nodes
% ================================================================
fprintf('[Step 5] Fourier partial sum at quadrature nodes ...\n');
f_j = fourier_partial_sum(f_hat, x_j, NF);
f_j = real(f_j);   % f is real; discard numerical noise in imaginary part

%% ================================================================
%  Step 6 — GLOF basis matrix  Phi(j,l+1) = x_j^nu * L_l^(alpha)(y_j)
% ================================================================
Phi = glof_basis_matrix(x_j, y_j, alpha, beta, lambda, M);

%% ================================================================
%  Step 7 — Normalization constants gamma_l
% ================================================================
gamma_l = glof_normalization(alpha, beta, M);

%% ================================================================
%  Step 8 — GLOF coefficients  c_tilde = G^{-1} Phi^T W f_j
% ================================================================
c_tilde = glof_coefficients(f_j, Phi, chi_j, gamma_l);

fprintf('  c_tilde norm (first 5): ');
fprintf('%.3e  ', c_tilde(1:min(5,end)));
fprintf('\n');

%% ================================================================
%  Step 9 — Evaluate reconstruction on test grid
% ================================================================
R     = real(glof_evaluate(c_tilde, z_eval, alpha, beta, lambda));
E_inf = max(abs(f_true - R));
E_rms = sqrt(mean((f_true - R).^2));
fprintf('\nMax error  : %.4e\n', E_inf);
fprintf('RMS error  : %.4e\n', E_rms);

%% ================================================================
%  Plots — reconstruction & error
% ================================================================
figure('Name','GLOF Reconstruction','Position',[100 100 800 600]);

subplot(2,1,1);
plot(z_eval, f_true, 'b-',  'LineWidth',2, 'DisplayName','Exact'); hold on;
plot(z_eval, R,      'r--', 'LineWidth',2, 'DisplayName','GLOF');
xlabel('x'); ylabel('f(x)');
title(sprintf('f(x)=x^{%.1f}(-\\log x)^{%d},   N_F=%d, M=%d, Q=%d', r, k, NF, M, Q));
legend('Location','northeast'); grid on;

subplot(2,1,2);
semilogy(z_eval, abs(f_true - R) + 1e-16, 'k-', 'LineWidth',1.5);
xlabel('x'); ylabel('|f - R_{N_F,M}f|');
title(sprintf('Pointwise error   (max = %.2e)', E_inf));
grid on;

%% ================================================================
%  GLOF coefficient tail — diagnostics
% ================================================================
figure('Name','GLOF Coefficient Tail','Position',[920 100 500 300]);
semilogy(0:M, abs(c_tilde) + 1e-20, 'bo-', 'MarkerSize',4);
xlabel('\ell'); ylabel('|c_\ell|');
title('GLOF coefficient magnitude'); grid on;

%% ================================================================
%  Convergence test — vary NF
% ================================================================
fprintf('\n=== Convergence Test ===\n');
fprintf('%-5s %-5s %-5s   E_inf\n', 'NF','M','Q');

NF_list = [8, 16, 32, 64];
E_list  = zeros(size(NF_list));

for ii = 1:length(NF_list)
    NF_i  = NF_list(ii);
    lam_i = floor(theta * NF_i);
    bet_i = lam_i + 2*r;
    M_i   = floor(gamma_coeff * NF_i);
    Q_i   = Q_factor * NF_i;

    % Fourier coefficients for this NF
    n_i  = (-NF_i : NF_i)';
    fh_i = zeros(2*NF_i+1, 1);
    for jj = 1:length(n_i)
        n = n_i(jj);
        re = @(x) x.^r .* (-log(x)).^k .* cos(2*pi*n*x);
        im = @(x) x.^r .* (-log(x)).^k .* sin(2*pi*n*x);
        fh_i(jj) = integral(re,0,1,'RelTol',1e-12,'AbsTol',1e-14,'Waypoints',[1e-10,1e-3]) ...
                 - 1i*integral(im,0,1,'RelTol',1e-12,'AbsTol',1e-14,'Waypoints',[1e-10,1e-3]);
    end

    [y_i, om_i] = laguerre_gauss_nodes(Q_i, alpha);
    x_i   = exp(-y_i / (bet_i + 1));
    chi_i = x_i.^(lam_i - bet_i) .* (bet_i+1)^(-(alpha+1)) .* om_i;
    fj_i  = real(fourier_partial_sum(fh_i, x_i, NF_i));
    Phi_i = glof_basis_matrix(x_i, y_i, alpha, bet_i, lam_i, M_i);
    gam_i = glof_normalization(alpha, bet_i, M_i);
    c_i   = glof_coefficients(fj_i, Phi_i, chi_i, gam_i);
    R_i   = real(glof_evaluate(c_i, z_eval, alpha, bet_i, lam_i));

    E_list(ii) = max(abs(f_true - R_i));
    fprintf('%-5d %-5d %-5d   %.4e\n', NF_i, M_i, Q_i, E_list(ii));
end

figure('Name','Convergence','Position',[100 720 500 300]);
loglog(NF_list, E_list, 'bo-', 'LineWidth',2, 'MarkerSize',8); hold on;
% Reference slope lines
loglog(NF_list, E_list(1)*(NF_list/NF_list(1)).^-2, 'k--', 'DisplayName','slope -2');
loglog(NF_list, E_list(1)*(NF_list/NF_list(1)).^-3, 'k:',  'DisplayName','slope -3');
xlabel('N_F'); ylabel('Max Error');
title('GLOF Convergence vs N_F'); legend('GLOF','slope -2','slope -3'); grid on;
