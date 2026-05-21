function gam = glof_norm_const(n, alpha, beta)
%GLOF_NORM_CONST  Normalization gamma_n^{(alpha,beta)} of the GLOF basis.
%
%   gam = GLOF_NORM_CONST(n, alpha, beta) returns
%       gamma_n^{(alpha,beta)}
%           = Gamma(n + alpha + 1) / [ (beta+1)^(alpha+1) * Gamma(n+1) ].
%
%   It is the squared L^2(w_{alpha,lambda}) norm of S_n^{(alpha,beta,lambda)}
%   where w_{alpha,lambda}(x) = (-log x)^alpha * x^lambda.
%
%   The implementation works in log-space to remain accurate for moderately
%   large n or alpha.
%
%   Input
%       n     : non-negative integer (scalar or array) -- the polynomial index
%       alpha : Laguerre parameter (alpha > -1)
%       beta  : GLOF parameter (beta > -1)
%
%   Output
%       gam : same size as n

    gam = exp(gammaln(n + alpha + 1) - gammaln(n + 1) ...
              - (alpha + 1) * log(beta + 1));
end
