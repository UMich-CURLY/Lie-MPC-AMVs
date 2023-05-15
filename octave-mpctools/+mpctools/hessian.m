function hes = hessian(func, var1, var2, output, name)
% `hess = hessian(func, [var1=1], [var2=var1], [output], [name])`
%
% Returns a CasADi function that gives the hessian of `func` with respect
% to the inputs `var1` and `var2`. `output` decides which output to find the
% hessian of if `func` has multiple outputs. By default, `output`, is all of
% the outputs of `func`.
%
% Note that `var1`, `var2` and `output` can be strings giving the input/output
% name or integers giving the (1-based) position of the input/output. `output`
% can also be a cell array, at which point `hess` will be a CasADi function with
% multiple outputs.
%
% CasADi does not currently support direct calculations of mixed hessians, so
% if `var1` and `var2` are not the same, then this function will fail with an
% error message from CasADi. If you do need a mixed hessian, you should simply
% call `mpctools.jacobian()` twice. Hopefully this limitation will be removed
% in a future CasADi release.
%
% `name` gives the name to use for the returned function object.
%
% Note that this function does not accept keyword arguments. Pass empty
% matrices as arguments to use the default values.
narginchk(1, 5);
if nargin() < 2 || isempty(var1)
    var1 = 1;
end
if nargin() < 3 || isempty(var2)
    var2 = var1;
end
if nargin() < 4 || isempty(output)
    output = num2cell(1:func.n_out());
end
if nargin() < 5 || isempty(name)
    name = ['hess_', func.name()];
end

% Check input and output variables.
if iscell(var1) || iscell(var2)
    error('var1 and var2 cannot be cell arrays!');
end
[~, var1] = mpctools.getargindex_(func, var1, 'in');
[~, var2] = mpctools.getargindex_(func, var2, 'in');

if ~iscell(output)
    output = {output};
end
[~, outnames] = mpctools.getargindex_(func, output, 'out');

% Use factory function.
outnames = cellfun(@(out) hessname(out, var1, var2), outnames, ...
                   'UniformOutput', false());
funcargs = mpctools.names2cell_(func.name_in());
if ~isequal(var1, var2)
    % Issue warning.
    warning(['CasADi does not support direct evaluation of mixed hessians.', ...
             ' See `help mpctools.hessian` for more information.']);
end
hes = func.factory(name, funcargs, outnames);

end%function

function s = hessname(out, v1, v2)
    s = sprintf('hess:%s:%s:%s', out, v1, v2);
end%function

