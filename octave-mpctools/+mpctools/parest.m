function solver = parest(varargin)
% `solver = parest(f, h, par, data, N, ...)`
%
% Returns ControlSolver object for solving parameter estimation problems.
%
% **Note**: This function is experimental, and its interface may change in the
% future.
%
% Inputs are as follows:
%
% * `f`: Casadi Function
%   
%   > Gives model for evolution of the system. By default, `f` is assumed to be
%   in discrete-time, but continuous-time `f` can be used with collocation
%   (see argument `N`).
%
% * `h`: Casadi Function
%
%   > Gives measurement function.
%
% * `par`: Cell array of strings
%
%   > Cell array giving names of parameters to estimate. Names cannot overlap
%   with the standard variables "x", "y", "v", etc. Sizes are inferred from
%   function arguments, which means each parameter to be estimated must appear
%   in some function (most likely "f" or "h").
%
% * `data`: Struct of arrays
%
%   > Holds all problem data that is fixed for the optimization problem, e.g.,
%   values of measurements $y$ and fixed inputs $u$. It must at least include
%   the field "y".
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
%   > Give the bounds and intial guess for the system variables. Each entry of
%   these structs should have the time dimension last (e.g., field "x" should
%   be of size `[N.x, N.t + 1]`, and "u" should be `[N.u, N.t]`). Any fields
%   that aren't given default to $+\infty$, $-\infty$, and $0$ respectively.
%
%   > Bounds and guesses for estimated parameters can be speficied in these
%   structs.
%
%   > When using collocation, if these structs contain "x" entries and not "xc"
%   entries, then values for "xc" will be inferred using linear interpolation.
%   If you do not want this behavior, you need to explicitly provide the "xc"
%   entry.
%
% * `Delta` : Scalar
%   
%   > Defines the timestep. Must be provided to use collocation.
%
% * `l` : Casadi Function
%
%   > Function giving the "stage cost" for the optimization problem. By default,
%   this is $l(v) = \|v\|^2$ with $v = h(x) - y$, i.e., the measurement error.
%   Note that this stage cost is applied at every time point, so it should not
%   include "u" as an argument.
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
% * `casaditype` : String [`'SX'`]
%
%   > Chooses which Casadi type to use for the NLP. If the functions are all
%   matrix-based or more complicated (e.g., Casadi integrators), or if the NLP
%   is very large then `'MX'` is usually best. If the functions are fairly
%   simple (e.g., a sequence of scalar operations) and if the NLP is small to
%   medium size, then `'SX'` is typically better. Default is `'SX'`.
%
% * `g` : Casadi Function
%
%   > If supplied, indicates that the model is a DAE system of the form
%   $x^+ = f(x,z,u)$, $g(x,z) = 0$ with differential states $x$ and algebraic
%   states $z$.
%
%   > As with any other function, you can define `g` to take any set of
%   arguments. Note that if `g` takes `u` as an argument, you should supply an
%   extra value of `u` (i.e., give `N.t + 1` values of `u`) to use for the final
%   constraint. Otherwise, the constraint at time `N.t` is written using
%   `u(N.t - 1)`.
%
% * `e` : Casadi Function
%
%   > Function that defines path constraints $e(...) \leq 0$. Note that the
%   arguments can be modified using `funcargs.e`. No constraint is written for
%   the terminal time point.
%
%   > To soften these constraints, include an "s" entry in `N` to give the
%   number of slacks you need. Note that if `funcargs.e` does not contain "s",
%   then the constraints are written as $e(x,u) \leq s$, which means `N.s` must
%   correspond to the number of components of `e`. Otherwise, you may use "s"
%   as an explicit argument to `e` (e.g., if you want to only soften certain
%   constraints).
%
% * `singleshooting` : Logical [`false`]
%
%   > Specifies whether to use single shooting to remove the system model from
%   the NLP, which trades problem size for sparsity. Note that single shooting
%   cannot be used with collocation.
%
% The output is a ControlSolver object.

% Handle arguments.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('f', 'required', 'casadifunc');
    parser.add('h', 'required', 'casadifunc');
    parser.add('par', 'required', {'cell', 'row'});
    parser.add('data', 'required', 'struct');
    parser.add('N', 'required', 'struct');
    parser.add('lb', struct(), 'struct');
    parser.add('ub', struct(), 'struct');
    parser.add('guess', struct(), 'struct');
    parser.add('Delta', [], {'scalar', 'pos'});
    parser.add('l', [], 'casadifunc');
    parser.add('funcargs', struct(), 'struct');
    parser.add('casaditype', 'SX', 'str');
    parser.add('g', [], 'casadifunc');
    parser.add('e', [], 'casadifunc');
    parser.add('singleshooting', false(), 'bool');
    parser.add('**kwargs');
end
[args, extraargs] = parser.parse(varargin{:});

% Create variables, bounds, and guesses.
N = args.N;
if ~all(isfield(N, {'x', 'y', 't'}))
    error('N must contain fields x, y, and t!');
elseif N.t <= 0
    error('N.t must be >= 1!');
end
usecollocation =  isfield(N, 'c') && N.c > 0;
if ~isfield(N, 'v')
    N.v = N.y;
end

% Set default objective function.
if isempty(args.l)
    args.l = mpctools.getCasadiFunc(@(v) v'*v, [N.v], {'v'}, {'l'});
end

% Figure out function arguments.
funcargs = mpctools.readfuncargs_({'f', 'h', 'l', 'g', 'e'}, args);
funcargs = mpctools.structupdate(funcargs, args.funcargs);

% Add slack to the constraint function if needed.
slackcon = mpctools.structget(N, 's', 0) > 0;
if slackcon && ~isempty(args.e) && ~ismember('s', funcargs.e)
    args.e = mpctools.addslackvar_(args.e, 's', -1, args.casaditype);
    funcargs.e = [funcargs.e, {'s'}];
end

% Wrap h if v is not an explicit variable.
if ~ismember(funcargs.h, 'v')
    args.h = mpctools.addslackvar_(args.h, 'v');
    funcargs.h = [funcargs.h, {'v'}];
end

% Check sizes for all function arguments.
[argsizes, allargs] = mpctools.readargsizes_(funcargs, args);

% Get sizes of all variables.
varlist = {'x', 'v'};
if usecollocation
    varlist = [varlist, {'xc'}];
end
if ~isempty(args.g)
    if ~isfield(N, 'z')
        error('Must provide N.z to use g!');
    end
    varlist = [varlist, {'z'}];
    if usecollocation
        varlist = [varlist, {'zc'}];
    end
end
if ~isfield(args.data, 'y')
    error('Must provide field "y" in data struct!');
elseif ~isequal(size(args.data.y), [N.y, N.t + 1])
    error('data.y must be size [N.y, N.t + 1] (expected [%s] but got [%s]!', ...
          mpctools.row2str([N.y, N.t + 1]), ...
          mpctools.row2str(size(args.data.y)));
end
if slackcon
    varlist = [varlist, {'s'}];
    if ~isfield(args.lb, 's')
        args.lb.s = 0; % Slacks are nonnegative by default.
    end
end
varlist = [varlist, intersect({'absx', 'absv', 'absd'}, allargs)];

% Make variable and parameter symbols. Note that "variable" refers to
% optimization variable, which includes the parameters to be estimated.
[sym, varlist, parlist] = mpctools.getocpsymbols(args.casaditype, argsizes, ...
                                                 N, varlist, args.par, ...
                                                 args.data);
if args.singleshooting
    sym = mpctools.singleshooting(sym, args.f, funcargs.f);
end

% Build struct of keyword arguments for controlproblem_.
kwargs = struct();
kwargs.N = N;
kwargs.funcs = struct('f', args.f, 'l', args.l, 'h', args.h, 'g', args.g, ...
                      'e', args.e);
kwargs.funcargs = funcargs;
kwargs.lb = args.lb;
kwargs.ub = args.ub;
kwargs.con = {};
kwargs.guess = args.guess;
kwargs.sym = sym;
kwargs.var = varlist;
kwargs.par = parlist;
kwargs.parval = args.data;
kwargs.Delta = args.Delta;
kwargs.obj = 0;
kwargs.singleshooting = args.singleshooting;
kwargs.finalstagecost = true();

% Call subfunction.
solver = mpctools.controlproblem_('**', kwargs, '**', extraargs);

end%function

