function integrator = getCasadiIntegrator(varargin)
% `integrator = getCasadiIntegrator(f, Delta, argsizes, argnames, ...)`
%
% Returns a Casadi `Integrator` object.
%
% Arguments are as follows:
%
% * `f` : Function
%
%   > Function that defines the ODE right-hand side.
%   Note that the first argument of `f` must be the differential variables $x$.
%
% * `Delta` : Positive Scalar
%
%   > Integration timestep.
%
% * `argsizes` : Cell or Vector
%
%   > Gives the size of each argument to `f`
%
% * `argnames` : Cell of Strings
%
%   > Gives names for each argument.
%
% * `funcname` : String
%
%   > Name for the Integrator object. Must be a valid variable identifier.
%
% * `wrap` : Logical [`true`]
%
%   > Whether to wrap the Integrator object so that it can be called via
%   $f(x, u, p, ...)$. If `false`, it must be called as
%   `f('x0', x, 'p', vertcat(u, p, ...))`, so you will typically always want
%   `wrap` to be `true`.
%
% * `Nt` : Integer [1]
%
%   > Number of time points to include. If `Nt` is larger than 1, all of the
%   arguments of `f` must be vectors. The returned function will then take
%   a vector for $x$ and matrices for all other parameters (with time along the
%   second dimension), and will return a matrix of the next `Nt` states.
%
%   > Note that if `wrap` is `false`, this argument has no effect.
%
% * `options` : Struct
%
%   > Struct of options to send to cvodes.
%
% * `solver` : String [`'cvodes'`]
%
%   > Which solver to use for integration. Typical options are `'cvodes'`,
%   `'idas'`, `'rk'`, and `'collocation'`. Consult the CasADi documentation for
%   more information about these options.
%
% * `casaditype` : String [`'SX'`]
%
%   > String `'SX'` or `'MX'` to decide which type of CasADi symbolic variables
%   to use in the ODE expression. Consult the CasADi User Guide for more
%   details.
%
%   > Note that previous versions allowed the use of a logical keyword argument
%   `scalar`, with `scalar=true` indicating `'SX'` and `scalar=false` indicating
%   `'MX'`. This option has been removed and should be replaced by `casaditype`.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('f', 'required', 'func');
    parser.add('Delta', 'required', {'scalar', 'pos'});
    parser.add('argsizes', 'required', 'row');
    parser.add('argnames', []);
    parser.add('funcname', 'int_f', {'str', 'varname'});
    parser.add('wrap', true(), 'bool');
    parser.add('Nt', 1, {'scalar', 'int', 'pos'});
    parser.add('options', struct(), 'struct');
    parser.add('solver', 'cvodes', 'str');
    parser.add('casaditype', 'SX', 'str');
end
args = parser.parse(varargin{:});
if ~casadi.has_integrator(args.solver)
    warning('Solver ''%s'' is unknown! Will lead to a CasADi error!', args.solver);
end

% Get symbolic expressions.
symbols = mpctools.getCasadiFunc_(args.f, args.argsizes, args.argnames, ...
                                  args.funcname, args.casaditype);
ode = struct();
ode.x = symbols.args{1}; % The differential variable.
ode.p = vertcat(symbols.args{2:end}); % Other variables (parameters).
ode.ode = symbols.fexpr; % Right-hand side.

intoptions = struct();
intoptions.verbose = false();
intoptions.inputs_check = false();
intoptions.disable_internal_warnings = true();
intoptions.regularity_check = false();
intoptions = mpctools.structupdate(intoptions, args.options);
intoptions.tf = args.Delta; % Save timestep to options struct.
integrator = casadi.integrator(args.funcname, args.solver, ode, intoptions);

% We have a raw integrator object, but often users want it to be in terms of
% their original variables.
if args.wrap
    if args.Nt > 1
        if any(cellfun(@(s) length(s) > 1 && s(2) > 1, symbols.sizes))
            error('All arguments must be vectors of Nt > 1!');
        end
        for i = 2:length(symbols.sizes)
            symbols.sizes{i} = [symbols.sizes{i}(1), args.Nt];
        end
    end
    inputvar = cell(1, length(symbols.names));
    for i = 1:length(inputvar)
        inputvar{i} = casadi.MX.sym(symbols.names{i}, symbols.sizes{i});
    end
    x0 = inputvar{1};
    if args.Nt == 1
        p = vertcat(inputvar{2:end});
        intout = integrator('x0', x0, 'p', p);
        xf = intout.xf;
        %TODO: allow matrix parameters.
    else
        outputvar = cell(args.Nt, 1);
        for t = 1:args.Nt
            p = cellfun(@(x) x(:,t), inputvar(2:end), 'UniformOutput', false());
            intout = integrator('x0', x0, 'p', vertcat(p{:}));
            x0 = intout.xf;
            outputvar{t} = x0;
        end
        xf = horzcat(outputvar{:});
    end
    integrator = casadi.Function(args.funcname, inputvar, {xf}, ...
                                 symbols.names, {args.funcname});
end

end%function

