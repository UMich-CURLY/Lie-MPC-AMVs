function x = cycle(x, xnew, dim)
% x = cycle(x, [xnew], [dim=2])
%
% Cycles array x along the given dimension. If xnew is given, it is used for the
% new slice. Otherwise, the first slice use used as the new one.
narginchk(1, 3);
if nargin() < 2
    xnew = [];
end
if nargin() < 3
    dim = 2;
end
if ndims(x) == 2 && dim == 2
    % Special case for 2D cycle along columns.
    if isempty(xnew)
        xnew = x(:,1);
    end
    x = [x(:,2:end), xnew];
else
    % General ND case.
    subs = repmat({':'}, 1, ndims(x));
    if isempty(xnew)
        subs{dim} = 1;
        xnew = subsref(x, substruct('()', subs));
    end
    subs{dim} = 2:size(x, dim);
    x = cat(dim, subsref(x, substruct('()', subs)), xnew);
end
end%function

