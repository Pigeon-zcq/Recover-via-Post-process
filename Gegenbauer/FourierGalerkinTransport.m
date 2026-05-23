function hf = FourierGalerkinTransport(N, p, q, Tfinal, cfl)
% FOURIERGALERKINTRANSPORT  Solve u_t - x*u_x = 0 by Fourier Galerkin + RK4.
%
% Problem (Chen & Shu 2015, Section 5, Example 6.5):
%   u_t - x*u_x = 0,  x in [-1,1],  t > 0
%   u(x,0) = (1+x)^(p/q)
%   periodic boundary conditions with period 2
%
% Method: Fourier Galerkin projection onto span{e^{ik*pi*x}, |k|<=N}.
% The ODE system for Fourier coefficients a_k(t) is:
%   d/dt a_k = sum_{j ~= k} [ j/(j-k) * (-1)^(j-k) ] * a_j
% (mass matrix from the operator -x*d/dx, see Fortran source).
%
% An exponential filter sigma_k = exp(-16*log(10)*(k/N)^r) is applied
% at each RK4 step to stabilise the computation.
%
% Inputs:
%   N      - Fourier truncation order (2N+1 modes)
%   p, q   - singularity exponent s = p/q
%   Tfinal - final time  (paper uses Tfinal = log(2))
%   cfl    - CFL number controlling dt  (paper uses cfl = 0.5)
%
% Output:
%   hf  - complex row vector of length 2N+1
%         hf(k+N+1) = a_k(Tfinal) for k = -N,...,N

    % Initial Fourier coefficients: a_k(0) = hat{u(x,0)}_k
    f_init = @(x) f_transport_ic(x, p, q);
    ak     = getPointValuefromFC(N, f_init);

    % Exponential filter
    a1    = 16 * log(10);
    r1    = floor(sqrt(2*N));
    k_idx = -N:N;
    sigma = exp(-a1 * (k_idx / N).^r1);

    % Mass matrix:  AM(i,j) = kj/(kj-ki)*(-1)^(kj-ki)  for ki ~= kj
    nmod = 2*N + 1;
    idx  = -N:N;
    AM   = zeros(nmod, nmod);
    for ii = 1:nmod
        ki = idx(ii);
        for jj = 1:nmod
            kj = idx(jj);
            if kj ~= ki
                AM(ii,jj) = kj / (kj - ki) * (-1)^(kj - ki);
            end
        end
    end

    % RK4 time integration  (ak is a row vector; transpose for matrix multiply)
    dt = cfl * Tfinal / (2*N + 1);
    t  = 0;
    while t < Tfinal
        if t + dt > Tfinal
            dt = Tfinal - t;
        end
        ak1 = (AM * ak.').' ;
        ak2 = (AM * (ak + 0.5*dt*ak1).').' ;
        ak3 = (AM * (ak + 0.5*dt*ak2).').' ;
        ak4 = (AM * (ak +     dt*ak3).').' ;
        ak  = ak + (dt/6) * (ak1 + 2*ak2 + 2*ak3 + ak4);
        t   = t + dt;
    end

    % Apply filter once at final time (not per-step, to avoid over-damping)
    hf = ak .* sigma;
end
