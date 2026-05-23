function y = transf_x2y(x_in, bet, a, b)
% TRANSF_X2Y  Coordinate transform: x in [a,b]  ->  y in [-1,1].
%
%   y = -1 + (x-a)^bet / ( 0.5*(b-a)^bet )
%
% Inverse of transf_y2x.  Corresponds to eq. (1.5) in Chen & Shu (2015).
%
% Inputs:
%   x_in - points in [a,b]  (array of any shape)
%   bet  - exponent = 1/q
%   a,b  - domain endpoints
%
% Output:
%   y    - corresponding points in [-1,1]  (same shape as x_in)

    mpsml      = 1e-16;
    half_b_bet = 0.5 * (b - a)^bet;
    y          = zeros(size(x_in));
    for i = 1:numel(x_in)
        if abs(x_in(i) - a) < mpsml
            y(i) = -1;
        elseif abs(b - x_in(i)) < mpsml
            y(i) = 1;
        else
            y(i) = -1 + (x_in(i) - a)^bet / half_b_bet;
        end
    end
end
