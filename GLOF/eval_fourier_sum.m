function val = eval_fourier_sum(x, ell, c)
%EVAL_FOURIER_SUM  Evaluate a truncated Fourier series.
%
%   val = EVAL_FOURIER_SUM(x, ell, c) returns the real-valued partial sum
%       sum_j  c(j) * exp(2*pi*i*ell(j)*x)
%   at every entry of x.  The imaginary part is discarded because the test
%   functions and their Fourier representations are real to roundoff.
%
%   Input
%       x   : array of evaluation points
%       ell : vector of mode indices, e.g. -N:N
%       c   : vector of Fourier coefficients (same length as ell)
%
%   Output
%       val : same size as x

    xrow = x(:).';
    ell  = ell(:);
    c    = c(:);

    % Compute the matrix exp(2*pi*i*ell*x) of size length(ell) x numel(x)
    expmat = exp((2*pi*1i) * ell * xrow);
    val    = real(c.' * expmat);
    val    = reshape(val, size(x));
end
