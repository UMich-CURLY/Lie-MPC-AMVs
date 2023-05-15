function [fcasadi, qcasadi] = getCasadiFunc(varargin)
% `[fcasadi, qcasadi] = getCasadiFunc(f, varsizes, [varnames], [funcname='f'], ...)`
%
% Returns a Casadi function using the function handle `f`.
%
% Inputs are as follows:
%
% * `f`: Function Handle
%
%   > Handle to function that is being transformed.
%
% * `varsizes`: Row Vector
%
%   > Gives the number of elements in each input argument (all must be vectors).
%
% * `varnames`: Cell Array
%
%   > List of strings for variable names. Defaults to generic names `x_1`,
%   `x_2`, etc.
%
% * `funcname`: String
%
%   > Name to use for function.
%
% * `rk4`: Logical [`false`]
%
%   > Whether the function should be discretized with an explicit RK4 method.
%   Note that the first argument of `f` is assumed to be the differential state.
%
% * `Delta`: Scalar
%
%   > Timestep to use if `rk4` is `true`.
%
% * `M`: Integer
%
%   > Number of iterations of RK4 to perform. Note that the total timestep is
%   always equal to `Delta`, so the stepsize for each iteration is `Delta/M`.
%
% * `quad` : Function Handle
%
%   > Handle to function that should be integrated over the interval. Note that
%   `quad` must accept the same arguments (in the same order) as `f`. If
%   provided, `rk4` must be `true`.
%
% * `quadname` : String
%
%   > Name to use for quadrature function. Default is `['Q', funcname]`.
%
% * `casaditype` : String [`'SX'`]
%
%   > String `'SX'` or `'MX'` to decide which type of CasADi symbolic variables
%   to use. In general `'SX'` should be used for primitive algebraic operations
%   (e.g., ODE right-hand sides), and `'MX'` should be used for higher-level
%   constructs. Consult the CasADi User Guide for more details.
%
%   > Note that previous versions allowed the use of a logical keyword argument
%   `scalar`, with `scalar=true` indicating `'SX'` and `scalar=false` indicating
%   `'MX'`. This option has been removed and should be replaced by `casaditype`.
%
% The outputs are casadi `Function` objects. Note that `qcasadi` is only
% defined if a `quad` argument was passed.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('f');
    parser.add('varsizes', 'required', {'row'});
    parser.add('varnames', []);
    parser.add('funcname', 'f', {'str', 'varname'});
    parser.add('rk4', false(), 'bool');
    parser.add('Delta', 1, {'scalar', 'pos'});
    parser.add('M', 1, {'scalar', 'pos'});
    parser.add('quad', []);
    parser.add('quadname', [], {'str', 'varname'});
    parser.add('casaditype', 'SX', 'str');
end
args = parser.parse(varargin{:});
needquad = ~isempty(args.quad);
if needquad && ~args.rk4
    error('quad requires rk4 to be true!');
elseif ~needquad && nargout() > 1
    error('Multiple outputs requires an input for quad!');
end

% Wrap the function if RK4 requested.
if args.rk4
    f = args.f;
    Delta = args.Delta;
    M = args.M;
    func = @(x, varargin) mpctools.rk4(args.f, x, varargin, args.Delta, args.M);
    if needquad
        quadfunc = @(x, varargin) mpctools.rk4(args.f, x, varargin, ...
                                               args.Delta, args.M, args.quad);
    end
else
    func = args.f;
end

% Get Casadi Function objects.
fcasadi = getfunc(func, args);
if needquad
    if isempty(args.quadname)
        args.quadname = ['Q', args.funcname];
    end
    args.funcname = args.quadname;
    qcasadi = getfunc(quadfunc, args);
else
    qcasadi = [];
end

end%function


function func = getfunc(func, args)
    % Gets function symbolics and returns Function object.
    
    % Call low-level function.
    symbols = mpctools.getCasadiFunc_(func, args.varsizes, args.varnames, ...
                                  args.funcname, args.casaditype);

    % Build Casadi function.
    func = casadi.Function(args.funcname, symbols.args, {symbols.fexpr}, ...
                           symbols.names, {args.funcname});
end%function

