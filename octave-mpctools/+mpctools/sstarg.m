function solver = sstarg(varargin)
% `solver = sstarg(f, h, l, N, lb, ub, guess, ...)`
%
% Returns ControlSolver object for solving MPC control problems.
%
% Inputs are as follows:
%
% * `f`: Casadi Function
%   
%   > Gives model for evolution of the system. By default, `f` is assumed to be
%   in discrete-time, but continuous-time `f` can be used with collocation
%   (see argument `N`).
%
% * `h` : Casadi Function
%
%   > Gives measurement function. If not given, the variable y is not included
%   in the optimization problem.
%
% * `l`: Casadi Function
%
%   > Gives the objective function for the system. If included, `l` must have
%   been defined with names, or funcargs.l must also be provided.
%
%   > If not given, a dummy objective is used so that any feasible solution is
%   also optimal.
%
% * `N`: Struct of Integers
%   
%   > Contains fields "x", "u", "y" and "t" to specify the dimension of $x$,
%   $u$, and $y$, as well as say how many time points to use. Optionally, it
%   can contain a "c" entry to say how many (interior) collocation points to use
%   (for a continuous-time model `f`).
%
% * `lb`: Struct
% * `ub`: Struct
% * `guess`: Struct
%   
%   > Give the bounds and initial guess for the system variables. Each entry of
%   these structs should have the time dimension last (e.g., field "x" should
%   be of size `[N.x, N.t + 1]`, and "u" should be `[N.u, N.t]`). Any fields
%   that aren't given default to $+\infty$, $-\infty$, and $0$ respectively.
%
%   > These structs can also contain "Du" fields to specify rate-of-change bounds
%   for the system.
%
% * `par` : Struct
%
%   > Defines values of fixed parameters. Special entries "xsp" and "usp" define
%   time-varying setpoint values (and must be sizes `[N.x, N.t + 1]` and
%   `[N.u, N.t]` respectively.
%
%   > Parameters can be vectors (including scalars) or matrices, and the values
%   can be time-varying constants. For time-varying vector (or scalar)
%   parameters, pass values as a matrix with each slice `p(:,t)` (i.e., each
%   column) giving the value at the corresponding time point. For time-varying
%   matrices, use a 3D array, with each slice `p(:,:,t)` giving the value at
%   each time. Note that access is modulo length, so if a parameter `p` is
%   passed as a 1 by 3 matrix, then its value will repeat every three time
%   points regardless of the horizon `N.t`.
%
% * `funcargs` : Struct
%
%   > Contains cell arrays for each function that define the sequence of
%   arguments for the given function.
%
%   > By default, function argument names are read directly from the
%   corresponding Casadi function. Thus, you only need to use `funcargs` if
%   you created the Casadi function without names, or if you used names that
%   are different from the names in the optimization problem.
%
% * `e` : Casadi Function
%
%   > Function that defines constraints for the system. If given, either
%   `e` must have been defined with names, or `funcargs.e` must be provided.
%
%   > To soften these constraints, include an "s" entry in `N` to give the
%   number of slacks you need. Note that if `funcargs.e` does not contain "s",
%   then the constraints are written as $e(x,u) \leq s$, which means `N.s` must
%   correspond to the number of components of `e`. Otherwise, you may use "s"
%   as an explicit argument to `e` (e.g., if you want to only soften certain
%   constraints).
%
% * `discretef` : Logical [`true`]
%
%   > If `true`, the model `f` is assumed to be in discrete time, and the
%   constraint is written as $f(x) = x$. If `false`, `f` is assumed to be
%   continuous time, and the constraint is $f(x) = 0$.
%
% * `udiscrete` : Logical Vector
%
%   > Vector of `true` and `false` to say whether components of $u$ should be
%   discretely (i.e., integer) valued. Note that to actually enforce this
%   restriction in the optimization problem, you will need to choose a solver
%   that supports discrete variables (e.g., bonmin).
%
% * `casaditype` : String [`'SX'`]
%
%   > Chooses which Casadi type to use for the NLP. Since `sstarg` problems are
%   typically quite small, `'SX'` is almost always the best choice. However, if
%   the model `f` includes any nonscalar operations (e.g., Casadi Integrator
%   calls), then you will need to us `'MX'`.
%
% * `g` : Casadi Function
%
%   > If supplied, indicates that the model is a DAE system of the form
%   `x^+ = f(x,z,u), g(x,z) = 0` with differential states `x` and algebraic
%   states `z`.
%
% * `customvar` : Cell array of strings
%
%   > List of custom variables to add to the problem. Note that each variable
%   must have been used in one of the functions so that sizes can be inferred.
%   Custom variables are also time-invariant, i.e., there is only one copy
%   that is used everywhere.
%
%   > Bounds on custom variables can be included in `lb` and `ub`. Otherwise,
%   variables are unbounded.
%
% The output is a ControlSolver object.

% Handle arguments.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('f', 'required', 'casadifunc');
    parser.add('h', [], 'casadifunc');
    parser.add('l', [], 'casadifunc');
    parser.add('N', 'required', 'struct');
    parser.add('lb', struct(), 'struct');
    parser.add('ub', struct(), 'struct');
    parser.add('guess', struct(), 'struct');
    parser.add('par', struct(), 'struct');
    parser.add('funcargs', struct(), 'struct');
    parser.add('e', [], 'casadifunc');
    parser.add('discretef', true(), 'bool');
    parser.add('udiscrete', [], {'logical', 'col'});
    parser.add('casaditype', 'SX', 'str');
    parser.add('g', [], 'casadifunc');
    parser.add('customvar', {}, {'cell', 'row'});
    parser.add('**kwargs');
end
[args, extraargs] = parser.parse(varargin{:});

% Create variables, bounds, and guesses.
N = args.N;
if ~isfield(N, 'x')
    error('N must contain field x!');
end
N.t = 1;

% Figure out which functions we need and what their arguments are.
funcs = struct();
funcs.f = args.f;
if ~isempty(args.l)
    funcs.l = args.l;
end
if ~isempty(args.h)
    if ~isfield(N, 'y')
        error('N.y must be given if l is given!');
    end
    funcs.h = args.h;
end
if ~isempty(args.e)
    funcs.e = args.e;
end
if ~isempty(args.g)
    if ~isfield(N, 'z')
        error('N.z must be given if g is given!');
    end
    funcs.g = args.g;
end

% Get default function arguments and then update from user.
funcargs = mpctools.readfuncargs_(fieldnames(funcs), funcs);
funcargs = mpctools.structupdate(funcargs, args.funcargs);

% Add slack to the constraint function if needed.
slackcon = mpctools.structget(N, 's', 0) > 0;
if slackcon && ~isempty(args.e) && ~ismember('s', funcargs.e)
    funcs.e = mpctools.addslackvar_(funcs.e, 's', -1, args.casaditype);
    funcargs.e = [funcargs.e, {'s'}];
end

% Parse function arguments.
[argsizes, allargs] = mpctools.readargsizes_(funcargs, funcs);

% Get sizes of everybody and make sure parameter shapes are correct.
varlist = {'x', 'u'};
if ~isempty(args.h)
    varlist = [varlist, {'y'}];
end
if ~isempty(args.g)
    varlist = [varlist, {'z'}];
end
absvars = intersect({'absx', 'absu', 'absy'}, allargs);
varlist = [varlist, absvars];
if slackcon
    varlist = [varlist, {'s'}];
end
parval = args.par;

% Make variable and parameter symbols. Note that we exclude the final x and y.
[sym, varlist, parlist] = mpctools.getocpsymbols(args.casaditype, argsizes, ...
                                                 N, varlist, args.customvar, ...
                                                 parval, false(), false());

% Write model constraints.
fargs = mpctools.getargs_(1, sym, funcargs.f);
steadystate = struct('lb', zeros(N.x, 1), 'ub', zeros(N.x, 1));
try
    steadystate.con = args.f(fargs{:});
catch err
    error('evaluating f(%s): %s', mpctools.row2str(funcargs.f), err.message);
end
if args.discretef
    steadystate.con = steadystate.con - sym.x{1};
end
steadystate.con = {steadystate.con};

% Build struct of keyword arguments for controlproblem_.
kwargs = struct();
kwargs.N = N;
kwargs.funcs = rmfield(funcs, 'f');
kwargs.funcargs = rmfield(funcargs, 'f');
kwargs.lb = args.lb;
kwargs.ub = args.ub;
kwargs.con = {'steadystate', steadystate};
kwargs.guess = args.guess;
kwargs.sym = sym;
kwargs.var = varlist;
kwargs.par = parlist;
kwargs.parval = parval;
kwargs.obj = 0;

% Adjust bounds for slackvars.
if ismember('s', varlist)
    kwargs.lb.s = zeros(N.s, 1);
    kwargs.ub.s = inf(N.s, 1);
end

% Add integer restriction to u.
if ~isempty(args.udiscrete)
    if length(args.udiscrete) ~= N.u
        error('udiscrete must be length N.u!');
    end
    kwargs.discreteness = struct('u', args.udiscrete);
end

% Call subfunction.
solver = mpctools.controlproblem_('**', kwargs, '**', extraargs);

end%function

