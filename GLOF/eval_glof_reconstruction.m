function R = eval_glof_reconstruction(x, g, alpha, beta, lambda)
%EVAL_GLOF_RECONSTRUCTION  Evaluate the GLOF partial sum.
%
%   R = EVAL_GLOF_RECONSTRUCTION(x, g, alpha, beta, lambda) returns
%       R(x) = sum_{n=0}^{m} g(n+1) * S_n^{(alpha,beta,lambda)}(x),
%   where m = length(g) - 1.  The result is real.
%
%   Input
%       x      : evaluation points (any shape), entries in [0,1]
%       g      : (m+1)-by-1 vector of GLOF coefficients
%       alpha  : Laguerre parameter
%       beta   : GLOF parameter
%       lambda : weight parameter
%
%   Output
%       R : same size as x, real-valued reconstruction.

    m = numel(g) - 1;
    S = glof_matrix(m, alpha, beta, lambda, x);     % (m+1) x numel(x)
    R = real(g(:).' * S);
    R = reshape(R, size(x));
end
