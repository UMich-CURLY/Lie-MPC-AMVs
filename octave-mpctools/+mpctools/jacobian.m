function jac = jacobian(func, var, output, name)
% `jac = jacobian(func, [var=1], [output], [name])`
%
% Returns a CasADi function that gives the derivative of `func` with respect
% to the input `var`. `output` decides which output to take the derivative of
% if `func` has multiple outputs.
%
% Note that `var` and `output` can both be strings giving the input/output name
% or integers giving the (1-based) position of the input/output. By default,
% `output` is all of the outpus of func if `var` is not a cell array, and 1
% if `var is a cell array.
%
% Either `var` or `output` can be cell arrays, at which point `jac` will be a
% CasADi function with multiple outputs. You can choose to either differentiate
% a single output with respect to multiple inputs, or multiple outputs with
% respect to a single input.
%
% `name` gives the name to use for the returned function object.
%
% Note that this function does not accept keyword arguments. Pass empty
% matrices as arguments to use the default values.
narginchk(1, 4);
if nargin() < 2 || isempty(var)
    var = 1;
end
if nargin() < 3 || isempty(output)
    if iscell(var)
        output = 1;
    else
        output = num2cell(1:func.n_out());
    end
end
if nargin() < 4 || isempty(name)
    name = ['jac_', func.name()];
end

% Check input and output variables.
if iscell(var) && iscell(output)
    error('Only 1 of var and output can be cell arrays!');
end
[~, innames] = mpctools.getargindex_(func, var, 'in');
[~, outnames] = mpctools.getargindex_(func, output, 'out');

% Use factory function.
if iscell(outnames)
    outnames = cellfun(@(out) jacname(out, innames), outnames, ...
                       'UniformOutput', false());
end
if ~iscell(innames)
    innames = {innames};
end
if ~iscell(outnames)
    outnames = cellfun(@(in) jacname(outnames, in), innames, ...
                       'UniformOutput', false());
end
funcargs = mpctools.names2cell_(func.name_in());
jac = func.factory(name, funcargs, outnames);

end%function

function s = jacname(out, in)
    s = sprintf('jac:%s:%s', out, in);
end%function

