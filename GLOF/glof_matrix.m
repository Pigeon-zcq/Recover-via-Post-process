function S = glof_matrix(m, alpha, beta, lambda, x)
%GLOF_MATRIX  Evaluate the GLOF basis S_n^{(alpha,beta,lambda)} at points x.
%
%   S = GLOF_MATRIX(m, alpha, beta, lambda, x) returns an (m+1)-by-numel(x)
%   matrix whose (n+1, q) entry is
%       S_n^{(alpha,beta,lambda)}(x(q))
%           = x(q)^{(beta-lambda)/2} * L_n^{(alpha)}( -(beta+1)*log x(q) ),
%   for n = 0,1,...,m.
%
%   Special case x = 0:
%       If theta := (beta - lambda)/2 > 0 the basis vanishes at x = 0.
%       This routine sets S_n(0) = 0 whenever theta > 0 and avoids
%       evaluating log(0).
%
%   Input
%       m      : highest GLOF index (m >= 0)
%       alpha  : Laguerre parameter, alpha > -1
%       beta   : GLOF parameter, beta > -1
%       lambda : weight parameter
%       x      : real vector of evaluation points, entries in [0,1]
%
%   Output
%       S : (m+1)-by-numel(x) matrix of basis values.

    theta = (beta - lambda) / 2;
    xrow  = x(:).';

    pos = (xrow > 0);

    % Argument y = -(beta+1)*log(x); use a safe dummy 0 where x <= 0
    % so that the Laguerre recurrence does not see Inf / NaN.
    y = zeros(size(xrow));
    y(pos) = -(beta + 1) * log(xrow(pos));

    L = laguerre_matrix(m, alpha, y);              % (m+1) x numel(x)

    xtheta = zeros(size(xrow));
    xtheta(pos) = xrow(pos).^theta;

    S = bsxfun(@times, xtheta, L);                 % (m+1) x numel(x)

    % At x = 0 the basis is zero if theta > 0 and equals L_n^{(alpha)}(inf)
    % which is not well-defined otherwise.  We use the convention
    % S_n(0) = 0 whenever theta > 0; if theta <= 0 (not the regime here)
    % we still return 0 by convention.
    S(:, ~pos) = 0;
end
