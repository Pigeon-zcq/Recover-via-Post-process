function L_mat = laguerre_poly_matrix(y, alpha, M)
% L_mat = laguerre_poly_matrix(y, alpha, M)
%
% Compute generalized Laguerre polynomial values via three-term recurrence.
%
% L_mat(j, l+1) = L_l^(alpha)(y(j)),   l = 0, 1, ..., M
%
% Recurrence (Section 3.2 of algorithm document):
%   L_0^(alpha)(y) = 1
%   L_1^(alpha)(y) = alpha + 1 - y
%   L_{l+1}^(alpha)(y) = [(2l+alpha+1-y)*L_l - (l+alpha)*L_{l-1}] / (l+1)
%
% Inputs:
%   y     - N x 1 column vector of evaluation points
%   alpha - Laguerre parameter
%   M     - highest polynomial degree
%
% Output:
%   L_mat - N x (M+1) matrix

y = y(:);
N = length(y);
L_mat = zeros(N, M + 1);

L_mat(:, 1) = 1.0;                        % L_0

if M >= 1
    L_mat(:, 2) = alpha + 1 - y;          % L_1
end

for l = 1 : M-1
    L_mat(:, l+2) = ((2*l + alpha + 1 - y) .* L_mat(:, l+1) ...
                     - (l + alpha)          .* L_mat(:, l  )) / (l + 1);
end
end
