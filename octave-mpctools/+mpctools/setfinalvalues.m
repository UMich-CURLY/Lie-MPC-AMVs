function x = setfinalvalues(x, t, y)
% x = setfinalvalues(x, t, y)
%
% Performs the assignment x(...,(end - t + 1):end) = y with the ellipsis being
% the appropriate number of colons (as would be very easy in numpy but not in 
% Matlab because in Matlab, N-dimensional arrays are second-class citizens).
%
% If y is the same size as x, then its final t values are used. If y is the
% appropriate size, then it is used directly. If y is a scalar, then it is
% used for every component.
narginchk(3, 3);

% Get size and indices of x.
xsize = size(x);
tkeep = xsize(end) - t;
start = prod(xsize(1:end-1))*tkeep + 1;

% Get size for y and grab appropriate values.
xsetsize = xsize;
xsetsize(end) = t;
ysize = size(y);
if isequal(ysize, xsize)
    y = y(start:end);
elseif isequal(ysize, xsetsize)
    y = y(:);
elseif ~isscalar(y)
    warning('y is the incorrect size. Assignment will fail!');
end

% Perform assignment.
x(start:end) = y;
end%function
