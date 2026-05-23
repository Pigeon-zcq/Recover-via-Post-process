function val = f_transport_ic(x, p, q)
% F_TRANSPORT_IC  Initial condition for u_t - x*u_x = 0.
%   u(x,0) = (1+x)^(p/q),  x in (-1,1)
    arg = x + 1;
    if arg >= 1e-16
        val = arg^(p/q);
    else
        val = 0;
    end
end
