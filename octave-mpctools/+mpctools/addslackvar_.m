function slackfunc = addslackvar_(func, name, mult, vartype)
% slackfunc = addslackvar_(func, [name='s'], [mult=1], [vartype])
%
% Adds a slack variable to the Casadi function func. mult specifies the
% multiplier to use for the new variable (typically +1 or -1), and vartype can
% be 'SX', or 'MX' to decide with Casadi type to use (if not specified, it is
% inferred).
%
% Only works for functions with a single output.
narginchk(2, 4);
if nargin() < 2
    name = 's';
end
if nargin() < 3
    mult = 1;
end
if nargin() < 4
    vartype = mpctools.getfunctiontype_(func);
end

% Make sure func only has a single output.
if func.n_out() ~= 1
    error('func can only have one output!');
end

% Get variables.
switch vartype
case 'SX'
    invar = func.sx_in();
    slack = @(s) casadi.SX.sym(name, s);
case 'MX'
    invar = func.mx_in();
    slack = @(s) casadi.MX.sym(name, s);
otherwise
    error('Unknown vartype "%s"', vartype);
end

% Build new function.
slackvar = slack(func.numel_out());
outvar = func(invar{:}) + mult*slackvar;
slackfunc = casadi.Function(func.name(), [invar, {slackvar}], {outvar});

end%function

