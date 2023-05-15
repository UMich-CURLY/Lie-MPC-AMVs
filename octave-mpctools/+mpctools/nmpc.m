function solver = nmpc(varargin)
% `solver = nmpc(f, l, N, x0, lb, ub, guess, ...)`
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
% * `l`: Casadi Function
%
%   > Gives the stage cost for the system. To include (possibly time-varying)
%   setpoints for $x$ and $u$, define `l` to take additional arguments called
%   "xsp" and "usp", and provide the values for these setpoints in the `par`
%   argument. To penalize rates of change, you may also use "Du" as an argument.
%
% * `N`: Struct of Integers
%   
%   > Contains fields "x", "u", and "t" to specify the dimension of $x$ and
%   $u$, as well as say how many time points to use. Optionally, it can contain
%   a "c" entry to say how many (interior) collocation points to use (for a
%   continuous-time model `f`).
%
% * `x0`: Vector  
%
%   > Gives the initial condition for the system. If not provided, no initial
%   condition is used.
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
%   > These structs can also contain "Du" and "Dx" fields to specify
%   rate-of-change bounds for the system.
%
%   > When using collocation, if these structs contain "x" entries and not "xc"
%   entries, then values for "xc" will be inferred using linear interpolation.
%   If you do not want this behavior, you need to explicitly provide the "xc"
%   entry.
%
% * `Vf`: Casadi Function
%   
%   > Gives the terminal cost for the system (zero if not given).
%
% * `Delta` : Scalar
%   
%   > Defines the timestep. Must be provided to use collocation.
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
% * `periodic` : Logical [`false`]
%
%   > Determines whether or not to add a periodicity constraint to the problem.
%
% * `e` : Casadi Function
%
%   > Function that defines path constraints $e(x,u) \leq 0$. Note that the
%   arguments can be modified using `funcargs.e`. One constraint is written for
%   each time point.
%
%   > To soften these constraints, include an "s" entry in `N` to give the
%   number of slacks you need. Note that if `funcargs.e` does not contain "s",
%   then the constraints are written as $e(x,u) \leq s$, which means `N.s` must
%   correspond to the number of components of `e`. Otherwise, you may use "s"
%   as an explicit argument to `e` (e.g., if you want to only soften certain
%   constraints).
%
% * `ef` : Casadi Function
%
%   > Defines the terminal constraint $e_f(x(N)) \leq 0$.
%
% * `udiscrete` : Logical Vector
%
%   > Vector of `true` and `false` to say whether components of $u$ should be
%   discretely (i.e., integer) valued. Note that to actually enforce this
%   restriction in the optimization problem, you will need to choose a solver
%   that supports discrete variables (e.g., bonmin).
%
% * `uprev` : Column Vector
%
%   > Gives the previous value of $u$ to calculate the first rate of change.
%
% * `discretel` : Logical [`true`]
%
%   > Indicates whether the objective function should be a discrete sum of
%   stage costs (`true`) or a collocation-based quadrature (`false`). Note that
%   the latter requires `N.c > 0`.
%
% * `casaditype` : String [`'SX'`]
%
%   > Chooses which Casadi type to use for the NLP. If the functions are all
%   matrix-based or more complicated (e.g., Casadi integrators), or if the NLP
%   is very large then `'MX'` is usually best. If the functions are fairly
%   simple (e.g., a sequence of scalar operations) and if the NLP is small to
%   medium size, then `'SX'` is typically better. Default is `'SX'`.
%
% * `xf` : Column Vector
%
%   > Gives the terminal value for `x`. If provided, a terminal equality
%   constraint is added. You may consider using a large terminal penalty instead
%   if you encounter solver issues.
%
% * `singleshooting` : Logical [`false`]
%
%   > Specifies whether to use single shooting to remove the system model from
%   the NLP, which trades problem size for sparsity. Note that single shooting
%   cannot be used with collocation.
%
% * `h` : Casadi Function
%
%   > If supplied, `y` variables are added to the formulation with the
%   constraint `y = h(x)`. `y` can then be used like any other variable (e.g.,
%   supplying bounds, using it in `l`, etc).
%
%   > Note that if `h` is specified, you also must specify `N.y` in `N`.
%
% * `finaly` : Logical [`true`]
%
%   > Decides whether to include `y` at the final time point or not. Default is
%   `true`, which means it is included and can be used in a terminal constraint.
%   Note that if `h` takes `u` as an argument, then the first `u` will be used,
%   and so unexpected behavior may result.
%
% * `g` : Casadi Function
%
%   > If supplied, indicates that the model is a DAE system of the form
%   `x^+ = f(x,z,u), g(x,z) = 0` with differential states `x` and algebraic
%   states `z`.
%
%   > As with any other function, you can define `g` to take any set of
%   arguments. Note that if `g` takes `u` as an argument, then the constraint
%   at time `N.t` is written using `u(N.t - 1)`.
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
    parser.add('l', 'required', 'casadifunc');
    parser.add('N', 'required', 'struct');
    parser.add('x0', [], {'numeric', 'col'});
    parser.add('lb', struct(), 'struct');
    parser.add('ub', struct(), 'struct');
    parser.add('guess', struct(), 'struct');
    parser.add('Vf', [], 'casadifunc');
    parser.add('Delta', [], {'scalar', 'pos'});
    parser.add('par', struct(), 'struct');
    parser.add('funcargs', struct(), 'struct');
    parser.add('periodic', false(), {'logical', 'scalar'});
    parser.add('e', [], 'casadifunc');
    parser.add('ef', [], 'casadifunc');
    parser.add('udiscrete', [], {'logical', 'col'});
    parser.add('uprev', [], {'numeric', 'col'});
    parser.add('discretel', true(), 'bool');
    parser.add('casaditype', 'SX', '/SX|MX/');
    parser.add('xf', [], {'numeric', 'col'});
    parser.add('singleshooting', false(), 'bool');
    parser.add('h', [], 'casadifunc');
    parser.add('finaly', true(), 'bool');
    parser.add('g', [], 'casadifunc');
    parser.add('customvar', {}, {'cell', 'row'});
    parser.add('**kwargs');
end
[args, extraargs] = parser.parse(varargin{:});

% Create variables, bounds, and guesses.
N = args.N;
if ~all(isfield(N, {'x', 'u', 't'}))
    error('N must contain fields x, u, and t!');
elseif N.t <= 0
    error('N.t must be >= 1!');
end
usecollocation =  isfield(N, 'c') && N.c > 0;

% Set default function arguments and then update from user input.
funcargs = mpctools.readfuncargs_({'f', 'l', 'Vf', 'e', 'ef', 'h', 'g'}, args);
funcargs = mpctools.structupdate(funcargs, args.funcargs);

% Add slack to the constraint function if needed.
slackcon = mpctools.structget(N, 's', 0) > 0;
if slackcon && ~isempty(args.e) && ~ismember('s', funcargs.e)
    args.e = mpctools.addslackvar_(args.e, 's', -1, args.casaditype);
    funcargs.e = [funcargs.e, {'s'}];
end

% Parse function arguments.
argsizes = mpctools.readargsizes_(funcargs, args);

% Get sizes of everybody and make sure parameter shapes are correct.
varlist = {'x', 'u'};
if usecollocation
    varlist = [varlist, {'xc'}];
end
if ~isempty(args.g)
    if ~isfield(N, 'z')
        error('Must supply N.z to use g!');
    end
    varlist = [varlist, {'z'}];
    if usecollocation
        varlist = [varlist, {'zc'}];
    end
end
if any(isfield(argsizes, {'Du', 'absDu'})) || isfield(args.lb, 'Du') ...
        || isfield(args.ub, 'Du')
    if isempty(args.uprev) && ~isfield(args.par, 'uprev')
        error('uprev must be given to use rate-of-change penalties/bounds!');
    end
    varlist = [varlist, {'Du'}];
end
if any(isfield(argsizes, {'Dx', 'absDx'})) || isfield(args.lb, 'Dx') ...
        || isfield(args.ub, 'Dx')
    varlist = [varlist, {'Dx'}];
end
if slackcon
    varlist = [varlist, {'s'}];
    if ~isfield(args.lb, 's')
        args.lb.s = 0; % Slacks are nonnegative by default.
    end
end
if ~isempty(args.h)
    if ~isfield(N, 'y')
        error('Must supply N.y to use h!');
    end
    varlist = [varlist, {'y'}];
end
absvars = {'absx', 'absu', 'absDu'};
absvars = absvars(isfield(argsizes, absvars));
if all(isfield(argsizes, {'xc', 'absx'})) && ~args.discretel
    % Need to include absxc.
    absvars = [absvars, {'absxc'}];
end
varlist = [varlist, absvars];
par = args.par;
if ~isempty(args.uprev)
    par.uprev = args.uprev;
end

% Get struct of Casadi symbols.
[sym, varlist, parlist] = mpctools.getocpsymbols(args.casaditype, argsizes, ...
                                                 N, varlist, args.customvar, ...
                                                 par, true(), args.finaly);
if args.singleshooting
    sym = mpctools.singleshooting(sym, args.f, funcargs.f);
end

% Build struct of keyword arguments for controlproblem_.
kwargs = struct();
kwargs.N = N;
kwargs.funcs = struct('f', args.f, 'l', args.l, 'e', args.e, 'h', args.h, ...
                      'g', args.g);
kwargs.funcargs = funcargs;
kwargs.lb = args.lb;
kwargs.ub = args.ub;
kwargs.con = {};
kwargs.guess = args.guess;
kwargs.sym = sym;
kwargs.var = varlist;
kwargs.par = parlist;
kwargs.parval = par;
kwargs.Delta = args.Delta;
kwargs.discretel = args.discretel;
kwargs.singleshooting = args.singleshooting;

% Add terminal cost and constraint.
if ~isempty(args.Vf)
    Vfargs = mpctools.getargs_(N.t + 1, sym, funcargs.Vf);
    try
        kwargs.obj = args.Vf(Vfargs{:});
    catch err
        error('evaluating Vf(%s): %s', mpctools.row2str(funcargs.Vf), ...
              err.message);
    end
end
if ~isempty(args.ef)
    efargs = mpctools.getargs_(N.t + 1, sym, funcargs.ef);
    termcon = struct();
    termcon.con = {args.ef(efargs{:})};
    termcon.lb = -inf(size(termcon.con{1}));
    termcon.ub = zeros(size(termcon.lb));
    kwargs.con = [kwargs.con, {'termcon', termcon}];
end

% Add periodic constraint.
if args.periodic
    periodic = struct();
    periodic.con = {sym.x{1} - sym.x{end}};
    periodic.lb = zeros(N.x, 1);
    periodic.ub = zeros(N.x, 1);
    kwargs.con = [kwargs.con, {'periodic', periodic}];
end

% Adjust bounds for slackvars.
if ismember(varlist, 's')
    kwargs.lb.s = zeros(N.s, N.t);
    kwargs.ub.s = inf(N.s, N.t);
end

% Add integer restriction to u.
if ~isempty(args.udiscrete)
    if length(args.udiscrete) ~= N.u
        error('udiscrete must be length N.u!');
    end
    kwargs.discreteness = struct('u', repmat(args.udiscrete, 1, N.t));
end

% Call subfunction.
solver = mpctools.controlproblem_('**', kwargs, '**', extraargs);

% Add in initial condition.
if ~isempty(args.x0)
    if length(args.x0) ~= N.x
        error('Invalid input for x0!');
    end
    solver.fixvar('x', 1, args.x0);
end
if ~isempty(args.xf)
    if length(args.xf) ~= N.x
        error('Invalid input for xf!');
    end
    solver.fixvar('x', N.t + 1, args.xf);
end

end%function

