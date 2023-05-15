function solver = nmhe(varargin)
% `solver = nmhe(f, h, u, y, l, N, lx, x0bar, lb, ub, guess, ...)`
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
%   > Gives measurement function.
%
% * `u` : Matrix
% * `y` : Matrix
%
%  > Give the known values for $u$ and $y$. Should be sized `[N.u, N.t]` and
%  `[N.y, N.t + 1]` respectively. Both can also be given in the `par` struct,
%  but values given as arguments override those in `par`.
%
%  > Note that `u` is only needed if the model has inputs, while `y` must
%  always be given (either as an argument or in `par`).
%
% * `l`: Casadi Function
%
%   > Gives the stage cost for the system error.
%
% * `N`: Struct of Integers
%   
%   > Contains fields "x", "u", "y" and "t" to specify the dimension of $x$,
%   $u$, and $y$, as well as say how many time points to use. Optionally, it
%   can contain a "c" entry to say how many (interior) collocation points to use
%   (for a continuous-time model `f`). It can also contain a "d" entry to give
%   the number of disturbance model states, which allows `d` to be used as an
%   argument in `f`, `h`, and `l` (the change in `d`, `Dd` can also be used
%   in `l`).
%
% * `lx`: Casadi Function
%
%   > Gives the arrival cost for x0 (zero if not given). Note that if
%   priorupdate (see below) is anything other than 'none', then this argument
%   cannot be specified, as the default quadratic prior is used.
%
% * `x0bar`: Vector  
%
%   > Gives the prior value to use for $\bar{x}_0$ in `lx`. Can also be given
%   in `par`.
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
%   > When using collocation, if these structs contain "x" entries and not "xc"
%   entries, then values for "xc" will be inferred using linear interpolation.
%   If you do not want this behavior, you need to explicitly provide the "xc"
%   entry.
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
% * `wadditive` : Logical [`false`]
%
%   > Decides whether w is an additive disturbance, i.e., $x^+ = f(x, u) + w$ or
%   is explicitly included in the model, i.e., $x^+ = f(x, u, w)$.
%
%   > When collocation is used, this choice means that $w$ is the difference
%   between the state at the right collocation *endpoint* and the state at the
%   next time point. In terms of the problem variables, interval `k`'s
%   collocation points are given by
%
%       [x(:,k), xc(:,:,k), x(:,k + 1) - w(:,k)]
%
%   and the value of $w$ is the instantaneous jump at the very end of the
%   interval.
%
% * `penalizevN` : Logical [`true`]
%
%   > Decides whether the final measurement error $v(N)$ should be penalized. If
%   true, an extra $l(0,v(N))$ term is added to the objective function (i.e.,
%   $l(w,v)$ with $w = 0$. Note that if you have custom arguments for $l$, only
%   $w$ is set to zero.
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
% * `singleshooting` : Logical [`false`]
%
%   > Specifies whether to use single shooting to remove the system model from
%   the NLP, which trades problem size for sparsity. Note that single shooting
%   cannot be used with collocation.
%
% * `priorupdate` : String [`'none'`]
%
%   > Specifies the prior update to use. Available options are  as follows
%   (with $k$ referring to the initial time of the MHE problem):
%
%   > * `'none'` : No automatic updates to prior parameters.
%   > * `'filtering'` : Updates prior parameters using an EKF step applied to
%       $\hat{x}(k - N | k - N)$.
%   > * `'smoothing'` : Updates prior parameters using $\hat{x}(k - N | k - 1)$.
%       Note that a correction is added to prevent "double-counting" of the
%       data $y(k - N)$ through $y(k)$.
%   > * `'hybrid'` : Identical to 'filtering' except that
%       $\hat{x}(k - N | k - 1)$ is used instead of $\hat{x}(k - N | k - N)$.
%
%   > For any choice besides `'none'`, the arrival cost is of the form
%   $\ell_x(x) = (x - \bar{x}_0)P^{-1}(x - \bar{x}_0)$ with parameter `Pinv`
%   giving the quadratic weight $P^{-1}$, and parameter `x0bar` giving the
%   the minimum value $\bar{x}_0$. This function is supplied automatically,
%   which means the `lx` argument should not be specified. Initial values for
%   `x0bar` and `Pinv` must both be specified in the `par` struct. These
%   parameters are then updated when `MHESolver.saveestimate()` is called.
%
%   > Note that prior updates are not supported for DAE models.
%
%   > For a linear system with uncorrelated $w$ and $v$, the `'filtering'` and
%   `'smoothing'` updates are equivalent to the (time-varying) Kalman Filter.
%   `'hybrid'` is not equivalent to the Kalman Filter but is faster and can give
%   better results for nonlinear systems.
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
% The output is an MHESolver object (subclass of ControlSolver with some extra
% methods specific to MHE problems).

% Handle arguments.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('f', 'required', 'casadifunc');
    parser.add('l', 'required', 'casadifunc');
    parser.add('h', 'required', 'casadifunc');
    parser.add('u', [], 'numeric');
    parser.add('y', [], 'numeric');
    parser.add('N', 'required', 'struct');
    parser.add('lx', [], 'casadifunc');
    parser.add('x0bar', [], {'numeric', 'col'});
    parser.add('lb', struct(), 'struct');
    parser.add('ub', struct(), 'struct');
    parser.add('guess', struct(), 'struct');
    parser.add('Delta', [], {'scalar', 'pos'});
    parser.add('par', struct(), 'struct');
    parser.add('funcargs', struct(), 'struct');
    parser.add('wadditive', false(), {'scalar', 'logical'});
    parser.add('penalizevN', true(), {'scalar', 'logical'});
    parser.add('casaditype', 'SX', '/SX|MX/');
    parser.add('g', [], 'casadifunc');
    parser.add('e', [], 'casadifunc');
    parser.add('singleshooting', false(), 'bool');
    parser.add('priorupdate', 'none', ...
               '/filter|filtering|smooth|smoothing|hybrid|none/');
    parser.add('customvar', {}, {'cell', 'row'});
    parser.add('**kwargs');
end
[args, extraargs] = parser.parse(varargin{:});

% Create variables, bounds, and guesses.
N = args.N;
if ~all(isfield(N, {'x', 'y', 't'}))
    error('N must contain fields x, y, and t!');
elseif N.t <= 0
    error('N.t must be >= 1! Consider mpctools.nlfilter() instead.');
end
usecollocation =  isfield(N, 'c') && N.c > 0;
if ~isfield(N, 'w')
    N.w = N.x;
end
if ~isfield(N, 'v')
    N.v = N.y;
end

% Figure out function arguments.
funcargs = mpctools.readfuncargs_({'f', 'l', 'h', 'lx', 'g', 'e'}, args);
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

% Wrap f if wadditive is true.
if args.wadditive
    if N.x ~= N.w
        error('w and x must have the same size when wadditive is true!');
    end
    if ismember('w', funcargs.f)
        error('w must not be an argument of f when wadditive is true!');
    end
    if ~usecollocation
        args.f = mpctools.addslackvar_(args.f, 'w');
        funcargs.f = [funcargs.f, {'w'}];
    end
end

% Check if a prior update is requested.
usingpriorupdate = true();
args.priorupdate = lower(args.priorupdate);
switch args.priorupdate
case 'none'
    usingpriorupdate = false();
case {'filter', 'filtering'}
    args.priorupdate = 'filtering';
case {'smooth', 'smoothing'}
    args.priorupdate = 'smoothing';
case 'hybrid'
    args.priorupdate = 'hybrid';
otherwise
    error('Unknown prior update "%s"!', args.priorupdate);
end

% Define arrival cost for prior update.
if usingpriorupdate
    if ~isempty(args.lx)
        error('Cannot specify custom lx is priorupdate is specified!');
    end
    funcargs.lx = {'x', 'x0bar', 'Pinv'};
    args.lx = mpctools.getCasadiFunc( ...
                @(x, x0bar, Pinv) (x - x0bar)'*Pinv*(x - x0bar), ...
                {N.x, N.x, [N.x, N.x]}, funcargs.lx, {'lx'});
    if any(~isfield(args.par, {'x0bar', 'Pinv'}));
        error('Must specify parameters x0bar and Pinv to use priorupdate!');
    end
end

% Check sizes for all function arguments.
[argsizes, allargs] = mpctools.readargsizes_(funcargs, args);

% Get sizes of everybody and make sure parameter shapes are correct.
varlist = {'x', 'w', 'v'};
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
if isfield(N, 'd')
    varlist = [varlist, {'d'}];
    if any(isfield(argsizes, {'Dd', 'absDd'}))
        varlist = [varlist, {'Dd'}];
        if isfield(argsizes, 'absDd')
            varlist = [varlist, {'absDd'}];
        end
    end
end
if slackcon
    varlist = [varlist, {'s'}];
    if ~isfield(args.lb, 's')
        args.lb.s = 0; % Slacks are nonnegative by default.
    end
end
varlist = [varlist, intersect({'absx', 'absw', 'absv', 'absd', 'absDs'}, ...
                              allargs)];

parval = args.par;
if ~isempty(args.u)
    parval.u = args.u;
end
if ~isempty(args.y)
    parval.y = args.y;
elseif ~isfield(parval, 'y')
    error('y must be provided as an argument or a field of par!');
end
if ~isempty(args.x0bar)
    parval.x0bar = args.x0bar;
end

% Make variable and parameter symbols.
[sym, varlist, parlist] = mpctools.getocpsymbols(args.casaditype, argsizes, ...
                                                 N, varlist, args.customvar, ...
                                                 parval);
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
kwargs.parval = parval;
kwargs.Delta = args.Delta;
kwargs.obj = 0;
kwargs.solvertype = 'mhe';
kwargs.singleshooting = args.singleshooting;
kwargs.priorupdate = args.priorupdate;

% Calculate objective for final measurement error v(N).
if args.penalizevN
    largs = mpctools.getargs_(N.t + 1, sym, funcargs.l);
    lzeros = cellfun(@(s) zeros(size(s)), largs, 'UniformOutput', false());
    lmask = ismember(funcargs.l, 'w');
    largs(lmask) = lzeros(lmask);
    try
        kwargs.obj = kwargs.obj + args.l(largs{:});
    catch err
        error('evaluating l(%s): %s', mpctools.row2str(funcargs.l), ...
              err.message);
    end
end

% Add arrival cost.
if ~isempty(args.lx)
    lxargs = mpctools.getargs_(1, sym, funcargs.lx);
    try
        kwargs.obj = kwargs.obj + args.lx(lxargs{:});
    catch err
        error('evaluating lx(%s): %s', mpctools.row2str(funcargs.lx), ...
              err.message);
    end
end

% Include information about model if a prior update is specified.
if usingpriorupdate
    if isfield(sym, 'z')
        error('Prior update is not supported for DAE systems!');
    end
    model = struct();
    funcnames = {'f', 'h', 'l'};
    for i = 1:length(funcnames)
        k = funcnames{i};
        model.(k) = struct('func', kwargs.funcs.(k), 'args', {funcargs.(k)});
    end
    
    % Create a discrete-time function for collocation models.
    if usecollocation
        model.f.func = getDiscretization(args.Delta, model.f.func, ...
                                         model.f.args, argsizes, ...
                                         args.casaditype);
        if args.wadditive
            model.f.func = mpctools.addslackvar_(model.f.func, 'w', 1, ...
                                                 args.casaditype);
            model.f.args = [model.f.args, {'w'}];
        end
    end
    
    % Save to kwargs.
    kwargs.mhemodel = model;
end

% Add w as x jump variables if wadditive and collocation are true.
if args.wadditive && usecollocation
    kwargs.collocjump = struct('x', {sym.w});
end

% Call subfunction.
solver = mpctools.controlproblem_('**', kwargs, '**', extraargs);

end%function

function F = getDiscretization(Delta, f, fargs, argsizes, casaditype)
    % Returns a discretized version of the given ODE. If casaditype is SX, then
    % RK4 (with M = 100) is used. If casaditype is MX, then CVODES is used.
    narginchk(5, 5);
    argsizes = cellfun(@(k) argsizes.(k), fargs, 'UniformOutput', false());
    switch casaditype
    case 'SX'
        F = mpctools.getCasadiFunc(f, argsizes, fargs, ...
                                   'rk4', true(), 'M', 100, 'Delta', Delta);
    case 'MX'
        F = mpctools.getCasadiIntegrator(f, Delta, argsizes, fargs);
    otherwise
        error('Invalid casaditype!');
    end
end%function
