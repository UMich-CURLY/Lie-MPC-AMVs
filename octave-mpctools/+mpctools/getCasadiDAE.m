function [dae, argorder] = getCasadiDAE(varargin)
% `[dae, argorder] = getCasadiDAE(Delta, f, [g], ...)`
%
% Arguments are as follows:
%
% * `Delta` : Positive Scalar
%
%   > Integration timestep.
%
% * `f` : Casadi Function
%
%   > Function that defines the ODE right-hand side.
%
% * `g` : Casadi Function
%
%   > Function that defines the algebraic constraints for the system. Note that
%   when there is overlap between arguments of `f` and `g`, the sizes must be
%   consistent.
%
%   > If `g` is not given, then the integrator is just a normal ODE.
%
% * `funcname` : String [`'dae'`]
%
%   > Name for the Integrator object. Must be a valid variable identifier.
%
% * `diffstate` : Integer or String [1]
% * `algstate` : Integer or String [2]
%
%   > Integer giving the position of or string giving the name of the variables
%   that define the differential variables $x$ and algebraic variables $z$.
%
% * `options` : Struct
%
%   > Struct of options to send to the solver.
%
% * `solver` : String [`'idas'`]
%
%   > Which solver to use for integration. Typical options are `'idas'` and
%   `'collocation'`. Consult the CasADi documentation for more information about
%   these options. Note that not all integrators support DAEs.
%
% The returned value `dae` is a CasADi Integrator object that returns the values
% of `x` and `z` after `Delta` time units. Arguments are passed in the order
% given in `argorder`, which is arguments of `f` in order followed by arguments
% unique to `g`.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('Delta', 'required', {'scalar', 'pos'});
    parser.add('f', 'required', 'casadifunc');
    parser.add('g', [], 'casadifunc');
    parser.add('funcname', 'dae', {'str', 'varname'});
    parser.add('diffstate', 1);
    parser.add('algstate', 2);
    parser.add('options', struct(), 'struct');
    parser.add('solver', 'idas', 'str');
end
args = parser.parse(varargin{:});
if ~casadi.has_integrator(args.solver);
    warning('Solver ''%s'' is unknown! Will lead to a CasADi error!', args.solver);
end
isdae = ~isempty(args.g);

% Get function arguments and sizes.
funcs = struct('f', args.f, 'g', args.g);
if isdae
    funcnames = {'f', 'g'};
else
    funcnames = {'f'};
end
funcargs = mpctools.readfuncargs_(funcnames, funcs);
argsizes = mpctools.readargsizes_(funcargs, funcs);

argorder = funcargs.f;
if isdae
    onlyg = ~ismember(funcargs.g, funcargs.f);
    argorder = [argorder, funcargs.g(onlyg)];
end

% Get symbols for everybody.
symbols = struct();
for i = 1:length(argorder)
    n = argorder{i};
    symbols.(n) = casadi.MX.sym(n, argsizes.(n));
end

% Choose differential and algebraic variables and get expressions.
ispar = true(size(argorder));
casadidae = struct();

[casadidae.x, xind] = getvar(args.diffstate, symbols, argorder);
ispar(xind) = false();
fargs = getargs(funcargs.f, symbols);
casadidae.ode = funcs.f(fargs{:});

if isdae
    [casadidae.z, zind] = getvar(args.algstate, symbols, argorder);
    ispar(zind) = false();
    gargs = getargs(funcargs.g, symbols);
    casadidae.alg = funcs.g(gargs{:});
end

pargs = cellfun(@(x) x(:), getargs(argorder(ispar), symbols), ...
                'UniformOutput', false());
casadidae.p = vertcat(pargs{:});

% Create Integrator object and wrap as Function.
args.options.tf = args.Delta;
integrator = casadi.integrator(args.funcname, args.solver, casadidae, ...
                               args.options);
intin = struct('x0', casadidae.x, 'p', casadidae.p);
if isdae
    intin.z0 = casadidae.z;
end
intout = integrator.call(intin);

symbols.xplus = intout.xf;
outputs = {'xplus'};
if isdae
    symbols.zplus = intout.zf;
    outputs = [outputs, {'zplus'}];
end
dae = casadi.Function(args.funcname, symbols, argorder, outputs);

end%function

function args = getargs(names, symbols)
    args = cellfun(@(s) symbols.(s), names, 'UniformOutput', false());
end%function

function [v, ind, name] = getvar(k, symbols, argorder)
    [ind, name] = mpctools.getargindex_(argorder, k, 'cell');
    v = symbols.(name);
end%function

