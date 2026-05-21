function p = bary_eval(xeval, xnodes, fnodes, w)
%BARY_EVAL  Barycentric interpolation evaluation.
%
%   p = BARY_EVAL(xeval, xnodes, fnodes, w) evaluates the polynomial
%   interpolant that takes the values fnodes(j) at xnodes(j), using the
%   barycentric weights w.  The output p has the same shape as xeval.
%
%   The "second" barycentric form is used:
%       p(x) = sum_j ( w_j/(x-x_j) * f_j ) / sum_j ( w_j/(x-x_j) ).
%   If x coincides (to a small tolerance) with a node, the corresponding
%   nodal value is returned directly.
%
%   Input
%       xeval  : array of evaluation points
%       xnodes : 1-by-(N+1) (or column) vector of node abscissae
%       fnodes : nodal values, same length as xnodes
%       w      : barycentric weights, same length as xnodes
%
%   Output
%       p : same size as xeval, real or complex according to fnodes.

    xnodes = xnodes(:).';
    fnodes = fnodes(:).';
    w      = w(:).';

    sz   = size(xeval);
    xrow = xeval(:).';

    p   = zeros(1, numel(xrow), 'like', fnodes(1));
    tol = 1e-14;

    for q = 1:numel(xrow)
        x = xrow(q);
        diffs = x - xnodes;
        hit = find(abs(diffs) < tol, 1);
        if ~isempty(hit)
            p(q) = fnodes(hit);
        else
            tmp = w ./ diffs;
            p(q) = sum(tmp .* fnodes) / sum(tmp);
        end
    end

    p = reshape(p, sz);
end
