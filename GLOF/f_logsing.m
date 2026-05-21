function val = f_logsing(x, r, k)
%F_LOGSING  Evaluate f(x) = x^r * (-log x)^k with f(0) := 0.
%
%   val = F_LOGSING(x, r, k) returns the value of the test function
%       f(x) = x^r * (-log(x))^k,   x in [0,1],   r > 0,   k >= 0 integer
%   at every entry of the array x.  The convention f(0) = 0 is enforced
%   so that the routine is safe to call at the left endpoint.
%
%   Input
%       x : real array, entries assumed in [0,1]
%       r : positive real exponent of x
%       k : non-negative integer power of (-log x)
%
%   Output
%       val : same size as x, equal to x.^r .* (-log(x)).^k where x>0
%             and equal to 0 where x == 0.

    val = zeros(size(x));
    idx = (x > 0);
    val(idx) = x(idx).^r .* (-log(x(idx))).^k;
end
