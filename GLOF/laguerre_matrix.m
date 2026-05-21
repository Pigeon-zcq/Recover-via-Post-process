function L = laguerre_matrix(m, alpha, y)
%LAGUERRE_MATRIX  Generalized Laguerre polynomials by recurrence.
%
%   L = LAGUERRE_MATRIX(m, alpha, y) returns an (m+1)-by-numel(y) matrix
%   whose (n+1, q) entry is L_n^{(alpha)}(y(q)) for n = 0,1,...,m.
%
%   The standard three-term recurrence is used:
%       L_0^{(alpha)}(y) = 1
%       L_1^{(alpha)}(y) = -y + alpha + 1
%       (n+1) L_{n+1}^{(alpha)}(y) = (2n + alpha + 1 - y) L_n^{(alpha)}(y)
%                                    - (n + alpha) L_{n-1}^{(alpha)}(y)
%
%   Input
%       m     : highest polynomial degree to evaluate (m >= 0)
%       alpha : Laguerre parameter, alpha > -1
%       y     : real vector of evaluation points (any shape)
%
%   Output
%       L : (m+1)-by-numel(y) matrix; row n+1 holds L_n^{(alpha)}.

    y = y(:).';
    Q = numel(y);
    L = zeros(m+1, Q);
    L(1, :) = 1;
    if m >= 1
        L(2, :) = -y + alpha + 1;
    end
    for n = 1:m-1
        L(n+2, :) = ((2*n + alpha + 1 - y) / (n + 1)) .* L(n+1, :) ...
                  - ((n + alpha) / (n + 1)) * L(n, :);
    end
end
