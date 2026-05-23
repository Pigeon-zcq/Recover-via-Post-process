% diag3_ex65.m — test varying nt and filter-per-step for Example 6.5

addpath('/Users/zcq/PaperCode/Recover/Gegenbauer');
clear; clc;

p = 1; q = 2; bet = 0.5;
a_sub = -0.5; b_sub = 0.0;
Tfinal = log(2);
f_et = @(x) sqrt(max(0,1+2*x));

fprintf('=== Effect of nt on reconstruction quality ===\n');
fprintf('N=10, lam=1, m=2, Galerkin hf\n');
N = 10; lam = 1; m = 2; nzn = N;
hf_gal = FourierGalerkinTransport(N, p, q, Tfinal, 0.5);
gt = compute_gt(lam, 1000);

fprintf('  %6s  %10s  hg(0:2)\n', 'nt', 'maxerr');
for nt = [30, 90, 300, 1000, 3000]
    hg = gcoe(hf_gal, N, m, lam, bet, a_sub, b_sub, gt, nt);
    [u, xx] = gegenbauer_reconstruct(hg, m, lam, bet, a_sub, b_sub, nzn, gt);
    ex = arrayfun(f_et, xx);
    fprintf('  %6d  %10.4e   [%.4f %.4f %.4f]\n', nt, max(abs(u-ex)), hg(1), hg(2), hg(3));
end

fprintf('\n=== Effect of filter-per-step on Galerkin hf ===\n');
fprintf('N=40, lam=1, m=10\n');
N=40; lam=1; m=10; nzn=N; nt=3*N;
gt = compute_gt(lam, max(m,nt)+1);

% Version 1: filter once at end (current state)
hf1 = FourierGalerkinTransport(N, p, q, Tfinal, 0.5);
hg1 = gcoe(hf1, N, m, lam, bet, a_sub, b_sub, gt, nt);
[u1, xx1] = gegenbauer_reconstruct(hg1, m, lam, bet, a_sub, b_sub, nzn, gt);
ex1 = arrayfun(f_et, xx1);
fprintf('Filter at end:   hg(0:3)=[%.4f %.4f %.4f %.4f], err=%.4e\n', ...
    hg1(1),hg1(2),hg1(3),hg1(4), max(abs(u1-ex1)));

% Version 2: no filter at all
hf2 = transport_no_filter(N, p, q, Tfinal, 0.5);
hg2 = gcoe(hf2, N, m, lam, bet, a_sub, b_sub, gt, nt);
[u2, xx2] = gegenbauer_reconstruct(hg2, m, lam, bet, a_sub, b_sub, nzn, gt);
ex2 = arrayfun(f_et, xx2);
fprintf('No filter:       hg(0:3)=[%.4f %.4f %.4f %.4f], err=%.4e\n', ...
    hg2(1),hg2(2),hg2(3),hg2(4), max(abs(u2-ex2)));

% Version 3: filter every step (Fortran style)
hf3 = transport_filter_every_step(N, p, q, Tfinal, 0.5);
hg3 = gcoe(hf3, N, m, lam, bet, a_sub, b_sub, gt, nt);
[u3, xx3] = gegenbauer_reconstruct(hg3, m, lam, bet, a_sub, b_sub, nzn, gt);
ex3 = arrayfun(f_et, xx3);
fprintf('Filter per step: hg(0:3)=[%.4f %.4f %.4f %.4f], err=%.4e\n', ...
    hg3(1),hg3(2),hg3(3),hg3(4), max(abs(u3-ex3)));

fprintf('\n=== hf spectrum comparison (k=0 to 10) ===\n');
fprintf('  k  |hf_end|   |hf_none|  |hf_step|\n');
for k=0:min(10,N)
    fprintf('  %2d  %.4e  %.4e  %.4e\n', k, ...
        abs(hf1(k+N+1)), abs(hf2(k+N+1)), abs(hf3(k+N+1)));
end

function hf = transport_no_filter(N, p, q, Tfinal, cfl)
    f_init = @(x) f_transport_ic(x, p, q);
    ak = getPointValuefromFC(N, f_init);
    nmod = 2*N+1;
    idx = -N:N;
    AM = zeros(nmod,nmod);
    for ii=1:nmod
        ki=idx(ii);
        for jj=1:nmod
            kj=idx(jj);
            if kj~=ki
                AM(ii,jj)=kj/(kj-ki)*(-1)^(kj-ki);
            end
        end
    end
    dt = cfl*Tfinal/(2*N+1);
    t = 0;
    while t < Tfinal
        if t+dt > Tfinal, dt=Tfinal-t; end
        ak1=(AM*ak.').'; ak2=(AM*(ak+0.5*dt*ak1).').';
        ak3=(AM*(ak+0.5*dt*ak2).').'; ak4=(AM*(ak+dt*ak3).').';
        ak=ak+(dt/6)*(ak1+2*ak2+2*ak3+ak4);
        t=t+dt;
    end
    hf = ak;
end

function hf = transport_filter_every_step(N, p, q, Tfinal, cfl)
    f_init = @(x) f_transport_ic(x, p, q);
    ak = getPointValuefromFC(N, f_init);
    a1=16*log(10); r1=floor(sqrt(2*N));
    sigma=exp(-a1*((-N:N)/N).^r1);
    nmod=2*N+1; idx=-N:N;
    AM=zeros(nmod,nmod);
    for ii=1:nmod
        ki=idx(ii);
        for jj=1:nmod
            kj=idx(jj);
            if kj~=ki, AM(ii,jj)=kj/(kj-ki)*(-1)^(kj-ki); end
        end
    end
    dt=cfl*Tfinal/(2*N+1);
    t=0;
    while t<Tfinal
        if t+dt>Tfinal, dt=Tfinal-t; end
        ak1=(AM*ak.').'; ak2=(AM*(ak+0.5*dt*ak1).').';
        ak3=(AM*(ak+0.5*dt*ak2).').'; ak4=(AM*(ak+dt*ak3).').';
        ak=ak+(dt/6)*(ak1+2*ak2+2*ak3+ak4);
        ak=ak.*sigma;  % filter every step
        t=t+dt;
    end
    hf=ak;
end
