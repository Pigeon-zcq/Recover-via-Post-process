function val = f_transport_exact(x, p, q)
% F_TRANSPORT_EXACT  Exact solution of u_t - x*u_x = 0 at T = log(2).
%   u(x, log2) = (1 + 2x)^(p/q)  on x in [-0.5, 0]
%   (singularity moves from x=-1 to x=-0.5 under the flow x(t)=x0*e^t)
    arg = 1 + 2*x;
    if arg >= 1e-16
        val = arg^(p/q);
    else
        val = 0;
    end
end
