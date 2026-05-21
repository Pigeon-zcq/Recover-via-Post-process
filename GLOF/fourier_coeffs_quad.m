function c = fourier_coeffs_quad(f_handle, ell, reltol, abstol)
%FOURIER_COEFFS_QUAD  High-accuracy Fourier coefficients by quadrature.
%
%   c = FOURIER_COEFFS_QUAD(f_handle, ell) returns Fourier coefficients
%       c_ell = integral_0^1 f(x) * exp(-2*pi*i*ell*x) dx,
%   for every mode in ell, computed with MATLAB's adaptive quadrature.
%   f_handle must accept array x values in [0,1].
%
%   c = FOURIER_COEFFS_QUAD(f_handle, ell, reltol, abstol) uses the
%   supplied integral tolerances.
%
%   Input
%       f_handle : function handle accepting array x in [0,1]
%       ell      : vector of Fourier mode indices, e.g. (-N:N).'
%       reltol   : optional relative tolerance
%       abstol   : optional absolute tolerance
%
%   Output
%       c : numel(ell)-by-1 complex vector of Fourier coefficients

    if nargin < 2
        error('fourier_coeffs_quad:NotEnoughInputs', ...
              'Provide a function handle and Fourier mode vector ell.');
    end

    if ~isa(f_handle, 'function_handle')
        error('fourier_coeffs_quad:InvalidFunction', ...
              'The first input must be a function handle.');
    end
    if nargin < 3 || isempty(reltol)
        reltol = 1e-12;
    end
    if nargin < 4 || isempty(abstol)
        abstol = 1e-14;
    end

    ell = ell(:);
    c   = zeros(numel(ell), 1);

    for j = 1:numel(ell)
        l = ell(j);
        integrand = @(x) f_handle(x) .* exp(-2*pi*1i*l*x);
        c(j) = integral(integrand, 0, 1, ...
                        'RelTol', reltol, 'AbsTol', abstol, ...
                        'ArrayValued', true);
    end
end
