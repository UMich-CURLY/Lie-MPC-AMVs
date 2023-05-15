function obj = sumstagecosts(N, l, var, funcargs, obj)
% obj = sumstagecosts(N, l, var, funcargs, [obj])
%
% Returns an SX object for the objective function.
narginchk(4, 5);
if nargin() < 5
    obj = 0;
end
for i = 1:N
    args = mpctools.getargs_(i, var, funcargs);
    try
        obj = obj + l(args{:});
    catch err
        error('evaluating l(%s): %s', mpctools.row2str(funcargs), err.message);
    end 
end
end%function
