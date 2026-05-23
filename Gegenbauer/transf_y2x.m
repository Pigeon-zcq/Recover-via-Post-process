function x = transf_y2x(y, bet, a, b)
% TRANSF_Y2X  Coordinate transform: y in [-1,1]  ->  x in [a,b].
%
%   x = ( 0.5*(b-a)^bet * (1+y) )^(1/bet) + a
%
% This map makes the singularity (x-a)^(1/q) analytic as a function of y,
% where bet = 1/q.  Corresponds to eq. (1.5) in Chen & Shu (2015).
%
% Inputs:
%   y    - points in [-1,1]  (array of any shape)
%   bet  - exponent = 1/q
%   a,b  - domain endpoints
%
% Output:
%   x    - corresponding points in [a,b]  (same shape as y)

    mpsml      = 1e-16;
    half_b_bet = 0.5 * (b - a)^bet;
    x          = zeros(size(y));
    for i = 1:numel(y)
        if abs(1 + y(i)) < mpsml
            x(i) = a;
        elseif abs(1 - y(i)) < mpsml
            x(i) = b;
        else
            x(i) = (half_b_bet * (1 + y(i)))^(1/bet) + a;
        end
    end
end
