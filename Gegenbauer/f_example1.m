function val = f_example1(x, p, q)
% F_EXAMPLE1  u(x) = cos(x) + sin(x)*(1+x)^(p/q)
% Examples:
%   p=1, q=2 → cos(x) + sin(x)*sqrt(1+x)    [Chen & Shu 2015, Ex. 6.1]
%   p=1, q=3 → cos(x) + sin(x)*(1+x)^(1/3)  [Chen & Shu 2015, Ex. 6.2]
    if 1 + x >= 1e-16
        val = cos(x) + sin(x) * (1 + x)^(p/q);
    else
        val = cos(x);
    end
end
