classdef ControlSolver < handle
% Class for holding an NLP solver for optimal control problems. Note that the
% user should not create instances of this class directly (i.e., using the
% constructor), but should instead use the problem-specific interfaces `nmpc`,
% `nmhe`, and `sstarg` (all of which return `ControlSolver` objects).
%
% Public attributes are as follows:
%
% * `par`
% * `lb`
% * `ub`
% * `guess`
% * `conlb`
% * `conub`
% * `discreteness`
%
%   > Structs that hold current values of parameters, bounds, etc. All
%   time-varying entries have time along the final dimension.
%
%   > Since values can be changed directly by users, care should be taken to
%   avoid changing the size of any fields, or else there will be odd error
%   messages on the next call to `solve()`.
%
%   > Entries in `discreteness` should all either be true or false to say
%   whether the particular variable is discrete-valued or not. Note that not all
%   solvers support discrete variables.
%
% * `var`
%
%   > Read-only struct that holds the most recent optimal solution. Before
%   `solve()` is called, this field is empty.
%
% * `stats`
%
%   > A struct that contains information about the most recent call to
%   `solve()`. Note that for solvers other than IPOPT, this struct may be
%   empty (due to the Casadi interface not populating these values).
%
% * `status`
%
%   > A string containing the solver's most recent return status. When using
%   IPOPT, success is indicated by `'Solve_Succeeded'`. For other solvers, it
%   may not have a meaningful value.
%
% * `obj`
%
%   > Objective value of most recent optimization.
%
% * `verbosity` : Integer Between 0 and 12 [0]
% * `timelimit` : Positive Scalar (in seconds) [60]
% * `maxiter` : Positive Integer [$\infty$]
% * `isQP` : Logical [`false`]
%
%   > Solver-generic options that can be changed using `set_*()` functions.
%   These are intended so that the user does not need to remember the solver-
%   specific names for common options. However, not all solvers may provide
%   access to all of these options. For example, due to CasADi limitations,
%   when using the solver Bonmin, setting `verbosity=0` does not hide all
%   output.
%
%   > Note that all of these values can be passed as keyword arguments to
%   `nmpc`, `nmhe`, and `sstarg` to avoid having to call the `set_*()` methods
%   after the object has been built.
%
%   > To change any other solver options, you will need to use the `init()`
%   function and the solver-specific option names (see also `getoptions()`).
%   Note that CasADi does not provide access to all solver-specific options.

properties (GetAccess=public, SetAccess=public)
    par = struct();
    conlb = struct();
    conub = struct();
    lb = struct();
    ub = struct();
    guess = struct();
    discreteness = struct();
end%properties

properties (GetAccess=public, SetAccess=protected)
    var = [];
    stats;
    status = 'N/A';
    obj = NaN();
    horizon = 0;
    maxhorizon = 0;
    verbosity = 5;
    isQP = false();
    name = 'solver';
    timelimit = 60;
    maxiter = [];
    nlp = [];
    sol = [];
    singleshooting = false();
    solveroptions = [];
    varsym = struct();
    parsym = struct();
end%properties

properties (GetAccess=protected, SetAccess=protected)
    optimizer;
    varsizes;
    solver;
    needupdate_ = false();
    collocpts = [];
    Delta = 1;
    map = struct(); % Contains maps from var structs to single long vector.
    shootingfunc;
    defaults = struct(); % Stores default fields.
    consym = struct();
    objsym = [];
end%properties

methods (Access=public)
    function self = ControlSolver(varargin)
        % self = ControlSolver(...)
        %
        % Creates the solver object.
        persistent parser
        if isempty(parser)
            parser = mpctools.ArgumentParser();
            parser.add('horizon', 'required', {'scalar', 'nonnegative', 'int'});
            parser.add('var');
            parser.add('varlb');
            parser.add('varub');
            parser.add('varguess');
            parser.add('obj');
            parser.add('con');
            parser.add('conlb');
            parser.add('conub');
            parser.add('discreteness');
            parser.add('verbosity', 0);
            parser.add('timelimit', 60);
            parser.add('solver', []);
            parser.add('isQP', false());
            parser.add('maxiter', [], {'scalar', 'pos', 'int'});
            parser.add('name', 'solver');
            parser.add('noinit', false(), 'bool');
            parser.add('collocpts', [], {'numeric', 'row'});
            parser.add('Delta', 1, {'numeric', 'scalar', 'pos'});
            parser.add('par', [], 'struct');
            parser.add('parval', [], 'struct');
            parser.add('singleshooting', false(), 'bool');
        end
        args = parser.parse(varargin{:});

        % Save bounds fields.
        self.lb = args.varlb;
        self.ub = args.varub;
        self.conlb = args.conlb;
        self.conub = args.conub;
        self.guess = args.varguess;
        self.discreteness = args.discreteness;

        % Save sizes.
        self.horizon = args.horizon;
        self.maxhorizon = args.horizon;
        self.varsizes = struct();
        varnames = fieldnames(args.varguess);
        for i = 1:length(varnames)
            vn = varnames{i};
            self.varsizes.(vn) = size(args.varguess.(vn));
        end
        self.collocpts = args.collocpts;
        self.Delta = args.Delta;

        % Save symbolic expressions.
        self.varsym = args.var;
        self.parsym = args.par;
        self.consym = args.con;
        self.objsym = args.obj;

        % Store options.
        self.set_isQP(args.isQP);
        self.set_name(args.name);
        self.set_verbosity(args.verbosity);
        self.set_timelimit(args.timelimit);
        self.set_maxiter(args.maxiter);
        self.singleshooting = args.singleshooting;
        self.set_solver(args.solver);

        % Build nlp.
        self.map.x = mpctools.VarLayout(self.varsym);
        self.map.g = mpctools.VarLayout(self.consym);
        if ~isequal(size(args.obj), [1, 1])
            error('The objective function is nonscalar (size [%s])!', ...
                  mpctools.row2str(size(args.obj)));
        end
        self.nlp = struct('f', self.objsym, ...
                          'g', self.map.g.symvec(self.consym), ...
                          'x', self.map.x.symvec(self.varsym));
        if self.singleshooting
            shooting = vertcat(self.consym.model{:});
            self.shootingfunc = casadi.Function('singleshooting', ...
                                                {self.nlp.x}, {shooting});
        end
        if (~isempty(args.parval) && ~isempty(fieldnames(args.parval))) ...
                || (~isempty(args.par) && ~isempty(fieldnames(args.par)))
            if isempty(args.parval) || isempty(args.par)
                error('Both par and parval must be given!');
            end
            self.par = args.parval;
            self.map.p = mpctools.VarLayout(self.parsym);
            self.nlp.p = self.map.p.symvec(self.parsym);
        else
            self.par = [];
        end

        % Initialize.
        if args.noinit
            self.needupdate_ = true(); % Set changed flag.
            self.optimizer = []; % Will be initialized later.
        else
            self.init();
        end
    end%function

    function copiedsolver = copy(self)
        % `solver = self.copy()`
        %
        % Returns a copy of the ControlSolver object in its current state.
        narginchk(1, 1);

        % Store some fields to the arguments list.
        args = struct();
        args.horizon = self.maxhorizon;
        args.var = self.varsym;
        args.varlb = self.lb;
        args.varub = self.ub;
        args.varguess = self.guess;
        args.obj = self.objsym;
        args.con = self.consym;
        args.conlb = self.conlb;
        args.conub = self.conub;
        args.discreteness = self.discreteness;
        args.verbosity = self.verbosity;
        args.timelimit = self.timelimit;
        args.solver = self.solver;
        args.isQP = self.isQP;
        args.maxiter = self.maxiter;
        args.name = self.name;
        args.collocpts = self.collocpts;
        args.Delta = self.Delta;
        args.parval = self.par;
        if ~isempty(self.par)
            args.par = self.parsym;
        end
        args.singleshooting = self.singleshooting;

        % Create new object.
        args.noinit = true(); % Will initialize later.
        copiedsolver = mpctools.ControlSolver('**', args);

        % Initialize using solver options and set a few other fields.
        copiedsolver.init(self.solveroptions);
        if self.horizon ~= self.maxhorizon
            copiedsolver.truncatehorizon(self.horizon);
        end
        if ~isempty(self.sol)
            copiedsolver.import_solution(self.sol.x, self.stats, self.sol);
        end
    end%function

    function init(self, varargin)
        % `self.init(options)`
        % `self.init('key1', value1, ['key2', value2], ...)`
        %
        % Initializes solver object using the provided `options` struct or
        % using multiple `'key', value` pairs. See Casadi user guide for more
        % information about available options.
        %
        % Note that calling `init()` saves the current state of the object as
        % the default state (i.e., the state that gets restored by calling
        % `reset()`).
        narginchk(1, inf());
        switch length(varargin)
        case 0
            options = struct();
        case 1
            options = varargin{1};
        otherwise
            options = struct(varargin{:});
        end

        % Update solver options.
        defopts = struct('solver', self.solver, 'isQP', self.isQP, ...
                         'verbosity', self.verbosity, ...
                         'maxiter', self.maxiter, ...
                         'timelimit', self.timelimit, ...
                         'discretevar', self.map.x.vec(self.discreteness));

        defopts = mpctools.genericsolveroptions('**', defopts);
        options = mpctools.structupdate(defopts, options);

        % Create solver.
        if casadi.has_nlpsol(self.solver)
            solverfunc = @(varargin) casadi.nlpsol(varargin{:});
        elseif casadi.has_conic(self.solver)
            if ~self.isQP
                error('Solver %s is a QP solver but self.isQP == false!', ...
                      self.solver);
            end
            solverfunc = @(varargin) casadi.qpsol(varargin{:});
        else
            error('Solver %s is not available!', self.solver);
        end
        self.solveroptions = options;
        self.optimizer = solverfunc(self.name, self.solver, self.nlp, options);

        % Save current values as defaults.
        self.savedefaults();

        % Shut off update flag.
        self.needupdate_ = false();
    end%function

    function solve(self)
        % `self.solve()`
        %
        % Solves the current optimization problem.

        % Update if parameters have changed.
        if self.needupdate_
            self.init();
        end

        % Get bounds, adjusting if the horizon is different.
        [varlb_, varub_] = self.getvarbounds();
        [conlb_, conub_] = self.getconbounds();

        % Solve NLP.
        args = struct('x0', self.map.x.vec(self.guess), ...
                      'lbx', self.map.x.vec(varlb_), ...
                      'ubx', self.map.x.vec(varub_),...
                      'lbg', self.map.g.vec(conlb_), ...
                      'ubg', self.map.g.vec(conub_));
        if ~isempty(self.par)
            args.p = self.map.p.vec(self.par);
        end
        args = mpctools.explodestruct(args);
        solution = self.optimizer(args{:});
        solutionstats = self.optimizer.stats();
        solutionvar = solution.x;
        self.import_solution(solutionvar, solutionstats, solution);
    end%function

    function import_solution(self, solvar, solstats, solsol)
        % `self.import_solution(var, [stats], [sol])`
        %
        % Imports a solution into the current object. This is useful, e.g., if
        % you want to solve the optimization problem using a solver outside of
        % CasADi.
        %
        % `var` should be either the struct solution (i.e., in the format of
        % `self.var`) or a single long vector of variables.
        %
        % `stats` should be struct giving solution statistics. In particular,
        % it should contain a "return_status" field. `sol` should be a struct
        % giving extra solution information, e.g., dual multipliers, etc., in
        % the form of `self.sol`. Note that both of these two arguments are
        % optional.
        narginchk(2, 4);
        if nargin() < 3
            solstats = struct();
        end
        if nargin() < 4
            solsol = struct();
        end

        % Save stats.
        self.stats = solstats;
        self.status = mpctools.structget(self.stats, 'return_status', '?');
        self.obj = full(mpctools.structget(solsol, 'f', NaN()));
        self.sol = solsol;

        % Reshape solution vector.
        if isstruct(solvar)
            solvec = [];
        else
            solvec = full(solvar);
            solvar = self.xvec2struct(solvec);
        end
        self.var = solvar;
        if self.singleshooting
            if isempty(solvec)
                solvec = self.map.x.ivec(self.var);
            end
            x = full(self.shootingfunc(solvec));
            self.var.x = [self.var.x, reshape(x, size(self.var.x, 1), [])];
        end

        % Add in NaNs for nonsense values if horizon has been truncated.
        if self.horizon < self.maxhorizon
            n = self.maxhorizon - self.horizon;
            varnames = fieldnames(self.var);
            for i = 1:length(varnames)
                v = varnames{i};
                self.var.(v) = mpctools.setfinalvalues(self.var.(v), n, NaN());
            end
        end

        % Add times.
        self.var.t = (0:self.horizon)*self.Delta;
        if ~isempty(self.collocpts)
            self.var.tc = bsxfun(@plus, self.var.t, ...
                                 self.Delta*self.collocpts(2:end-1)');
        end
    end%function

    function saveguess(self, newguess, toffset)
        % `self.saveguess([newguess], [toffset])`
        %
        % Stores a guess from a given struct or from the current solution.
        %
        % If `newguess` is provided, all of its fields are copied to the current
        % guess (`toffset` defaults to 0). If not, the current optimal solution
        % is used (`toffset` defaults to 1). Note that any extra fields in the
        % given guess are silently ignored. Also, any entries of `newguess` that
        % are NaN will not be used.
        %
        % `toffset` can be manually specified to explicitly say give the time
        % offset to use. Note that only the overlapping time points are actually
        % changed.
        narginchk(1, 3);
        usenew = (nargin() >= 2) && ~isempty(newguess);
        if ~usenew
            newguess = self.var;
        end
        if nargin() < 3 || isempty(toffset)
            if usenew
                toffset = 0;
            else
                toffset = 1;
            end
        end
        oldguess = self.guess; % Save this guy.

        % Now actually cycle everybody. Recall that the time dimension is always
        % last.
        fields = intersect(fieldnames(newguess), fieldnames(self.guess));
        if self.singleshooting
            fields = setdiff(fields, {'x'});
        end
        for i = 1:length(fields)
            f = fields{i};

            % Check sizes.
            guesssize = size(self.guess.(f));
            newsize = size(newguess.(f));
            if ~isequal(guesssize, newsize)
                error('Size mismatch for variable "%s"', f);
            end

            % Compute the time values.
            tmin = max(toffset, 0) + 1;
            tmax = min(guesssize(end), newsize(end) + toffset);

            tguess = (tmin - toffset):(tmax - toffset);
            tnew = tmin:tmax;

            % Actually perform the assignment.
            if ismatrix(self.guess.(f))
                % Matrix case is easy.
                self.guess.(f)(:,tguess) = newguess.(f)(:,tnew);
            else
                % ND case. A bit trickier.
                colons = repmat({':'}, 1, ndims(newguess.(f)) - 1);
                idxguess = substruct('()', [colons, {tguess}]);
                idxnew = substruct('()', [colons, {tnew}]);

                thisnew = subsref(newguess.(f), idxnew);
                self.guess.(f) = subsasgn(self.guess.(f), idxguess, thisnew);
            end
            
            % Screen out any NaNs. TODO: refactor to avoid this.
            naninds = isnan(self.guess.(f));
            if any(naninds(:))
                self.guess.(f)(naninds) = oldguess.(f)(naninds);
            end
        end
    end%function

    function fixvar(self, var, t, val, inds)
        % `self.fixvar(var, t, val, [inds])`
        %
        % Sets guess, lb, and ub for variable `var` at time `t` to `val`.
        %
        % If the variable is a vector, you can use `inds` to set only some of
        % the components.
        narginchk(4, 5);
        if nargin() < 5
            inds = ':';
        end
        if ~isscalar(t)
            error('t must be a scalar!');
        end
        fields = {'lb', 'ub', 'guess'};
        if isvector(val)
            idx = substruct('()', {inds, t});
        else
            idx = substruct('()', [repmat({':'}, 1, ndims(val)), {t}]);
        end
        for i = 1:length(fields)
            f = fields{i};
            self.(f).(var) = subsasgn(self.(f).(var), idx, val);
        end
    end%function

    function truncatehorizon(self, newhorizon)
        % `self.truncatehorizon(newhorizon)`
        %
        % Truncates the horizon of the optimization problem by removing
        % constraints and fixing variables to zero.
        %
        % Note that if there are terminal costs or constraints, these are not
        % shifted. Thus, if you really need the true shorter horizon problem,
        % you will need to call the original function again with the appropriate
        % horizon from the beginning.
        if ~isscalar(newhorizon) || round(newhorizon) ~= newhorizon
            error('newhorizon must be a scalar integer!');
        elseif newhorizon < 0
            error('newhorizon must be nonnegative!');
        elseif newhorizon > self.maxhorizon
            error('newhorizon must be <= the original horizon.');
        end
        self.horizon = newhorizon;
    end%function

    function problem = getQP(self, qpsolver, point)
        % `problem = getQP(self, [solver='gurobi'], [point])`
        %
        % Returns a struct of standard-form parameters to solve a quadratic
        % approximation of the current NLP.
        %
        % The `solver` argument chooses the output format. The default is
        % `'gurobi'`, which returns a struct with Gurobi-compatible fieldnames.
        % Also available are `'quadprog'`, which returns a struct for use with
        % Matlab's `quadprog`, and ``qp``, which returns a *cell array* of
        % arguments that can be given to Octave's `qp` function. Note that only
        % Gurobi supports discrete decision variables.
        %
        % If given, `point` is a struct that defines the point to use for
        % linearization. If not supplied, `self.guess` is used.
        %
        % Note that the objective function is the nominal objective function,
        % not the Lagrangian.
        narginchk(1, 3);
        if nargin() < 2
            qpsolver = [];
            point = [];
        elseif isstruct(qpsolver) && nargin() < 3
            point = qpsolver;
            qpsolver = [];
        elseif nargin() < 3
            point = [];
        end
        if isempty(qpsolver)
            qpsolver = 'gurobi';
        end
        if isempty(point)
            point = self.guess;
        end
        
        % Decision variable and parameters.
        x0 = self.map.x.vec(point);
        if isfield(self.nlp, 'p')
            p0 = self.map.p.vec(self.par);
        else
            p0 = [];
        end
        
        % Create an "oracle" function [x, p] -> [f, g].
        nlpfunc = casadi.Function('nlpfunc', self.nlp, {'x', 'p'}, {'f', 'g'});
        
        % Create a function for calculating derivative information. Note that
        % factory only returns the upper-triangular part of the hessian, so we
        % have to fill in the lower-triangular part manually.
        qpfunc = nlpfunc.factory('qpfunc', {'x', 'p'}, ...
                                 {'f', 'g', 'grad:f:x', 'hess:f:x:x', 'jac:g:x'});
        [f0, g0, grad_f0, hess_f0, jac_g0] = qpfunc(x0, p0);
        
        % Make matrices.
        problem = struct();
        
        problem.Q = sparse(hess_f0);
        problem.Q = problem.Q + triu(problem.Q, 1)';
        
        problem.obj = full(grad_f0) - problem.Q*x0;
        problem.objcon = -0.5*x0'*problem.Q*x0 - problem.obj'*x0 + full(f0);
        
        A = sparse(jac_g0);
        b = full(g0) - A*x0;
        [conlb_, conub_] = self.getconbounds();
        blb = self.map.g.vec(conlb_) - b;
        bub = self.map.g.vec(conub_) - b;
        eqcon = abs(blb - bub) < 1e-8;
        lbcon = ~eqcon & (blb ~= -inf());
        ubcon = ~eqcon & (bub ~= inf());
        ineqcon = lbcon | ubcon;
        beq = 0.5*(blb(eqcon) + bub(eqcon));

        % Get bounds and variable types.
        problem.lb = self.map.x.vec(self.lb);
        problem.ub = self.map.x.vec(self.ub);

        problem.vtype = repmat('C', size(A, 2), 1);
        problem.vtype(logical(self.map.x.vec(self.discreteness))) = 'I'; % Integers.
        
        % Tweak fields based on chosen solver.
        switch qpsolver
        case 'gurobi'
            problem.Q = problem.Q/2;
            problem.A = [A(eqcon,:); A(lbcon,:); A(ubcon,:)];
            problem.rhs = [beq; blb(lbcon); bub(ubcon)];
            problem.sense = [repmat('=', nnz(eqcon), 1); ...
                             repmat('>', nnz(lbcon), 1);
                             repmat('<', nnz(ubcon), 1)];
        case 'quadprog'
            intcon = find(problem.vtype == 'I');
            problem = struct('H', problem.Q, 'f', problem.obj, ...
                             'Aineq', [-A(lbcon,:); A(ubcon,:)], ...
                             'bineq', [-blb(lbcon); bub(ubcon)], ...
                             'Aeq', A(eqcon,:), 'beq', beq, ...
                             'lb', problem.lb, 'ub', problem.ub, ...
                             'options', optimoptions('quadprog'), ...
                             'solver', 'quadprog');
            if ~isempty(intcon)
                problem.intcon = intcon;
            end
        case 'qp'
            problem = {x0, problem.Q, problem.obj, A(eqcon,:), beq, ...
                       problem.lb, problem.ub, blb(ineqcon), ...
                       A(ineqcon,:), bub(ineqcon)};
        otherwise
            error('Unknown choice for solver!');
        end
    end%function

    function docstring = getoptions(self)
        % `self.getoptions()`
        %
        % Returns the options available for the current solver as a string.
        if casadi.has_nlpsol(self.solver)
            docstring = casadi.doc_nlpsol(self.solver);
        elseif casadi.has_conic(self.solver)
            docstring = casadi.doc_conic(self.solver);
        else
            error('Unknown solver: %s', self.solver);
        end
    end%function

    function s = xvec2struct(self, v)
        % `s = self.xvec2struct(v)`
        %
        % Reshapes the long vector `v` into a struct of named variables `s`.
        s = self.map.x.ivec(v);
    end%function

    function s = gvec2struct(self, v)
        % `s = self.gvectostruct(v)`
        %
        % Reshapes the long vector `v` into a struct of named constraints `s`.
        % Useful for multipliers on constraints.
        s = self.map.g.ivec(v);
    end%function

    function set_verbosity(self, val)
        maxverb = mpctools.MAX_VERBOSITY();
        self.verbosity = min([12, maxverb, max(0, round(val))]);
        self.needupdate_ = true();
    end%function

    function set_timelimit(self, val)
        if ~isscalar(val) || val <= 0
            error('Timelimit must be a positive scalar!');
        end
        self.timelimit = val;
        self.needupdate_ = true();
    end%function

    function set_name(self, val)
        if ~isvarname(val)
            error('Name "%s" is not a valid identifier!', val);
        end
        self.name = val;
        self.needupdate_ = true();
    end%function

    function set_isQP(self, val)
        if ~isscalar(val) || ~islogical(val)
            error('isQP must be a logical scalar!');
        end
        self.isQP = val;
        self.needupdate_ = true();
    end%function

    function set_maxiter(self, val)
        if ~isscalar(val) && ~isempty(val)
            error('maxiter must be a scalar!');
        end
        self.maxiter = max(1, round(val));
        self.needupdate_ = true();
    end%function

    function set_solver(self, solvername)
        % `self.set_solver([solvername])`
        %
        % Sets the solver to use for optimization.
        %
        % Note that setting the solver will reset any solver-specific options
        % that have been set via `init()`.
        narginchk(1, 2);
        if nargin() < 2
            solvername = [];
        end

        % Choose solver.
        if isempty(solvername)
            if self.isQP && casadi.has_conic('qpoases')
                solvername = 'qpoases';
            else
                solvername = 'ipopt';
            end
        elseif casadi.has_conic(solvername)
            if ~self.isQP
                warning('Solver %s is a QP solver but self.isQP == false!', ...
                        solvername);
            end
        elseif casadi.has_nlpsol(solvername)
            % Nothing to check.
        else
            warning(['Solver %s is not available! ', ...
                     'Will lead to error on init() or solve()!'], self.solver);
        end
        self.solver = solvername;
        self.solveroptions = struct();
        self.needupdate_ = true();
    end%function

    function reset(self)
        % Resets the object back to its default state.
        if ~all(isfield(self.defaults, self.defaultfields()))
            error('Default have not been saved! Forgot to self.init()?');
        end
        self.restoredefaults();
        self.init(self.solveroptions);
    end%function
    
    function cyclepar(self, varargin)
        % `self.cyclepar('par1', new1, ['par2, 'new2', ...])`
        %
        % Cycles one or more time-varying parameters. For each parameter given,
        % the oldest value is removed, everything is shifted forward by one
        % spot, and the new value is inserted.
        %
        % E.g., for parameter xsp, the call
        %```
        %    self.cyclepar('xsp', xspnew)
        %```
        % is equivalent to
        %```
        %    self.par.xsp = [self.par.xsp(:,2:end), xspnew]
        %```
        % This function works for time-varying scalars, vectors, and matrices.
        % Note that it assumes that the horizon is full (i.e.,
        % `self.horizon == self.maxhorizon`).
        narginchk(1, inf());
        if mod(length(varargin), 2) ~= 0
            error('Arguments must come in pairs!');
        end
        for i = 1:2:length(varargin)
            k = varargin{i};
            if ~ischar(k)
                error('Parameter names must be strings!');
            elseif ~isfield(self.par, k)
                error('Unknown parameter ''%s''!', k);
            end
            p = self.par.(k);
            if iscolumn(p)
                dim = 1;
            else
                dim = ndims(p);
            end
            newp = varargin{i + 1};
            self.par.(k) = mpctools.cycle(p, newp, dim);
        end
    end%function
end%methods

methods (Access=protected)
    function savedefaults(self)
        % Saves the current values of the object as the defaults.
        fields = self.defaultfields();
        self.defaults = struct();
        for i = 1:length(fields)
            f = fields{i};
            self.defaults.(f) = self.(f);
        end
    end%function

    function restoredefaults(self)
        % Restores the saved default values. Note that you need to call `init()`
        % to also reset any custom solver options.
        fields = self.defaultfields();
        for i = 1:length(fields)
            f = fields{i};
            self.(f) = self.defaults.(f);
        end
    end%function

    function fields = defaultfields(self)
        % Returns a list of the fields that are stored in defaults.
        fields = {'lb', 'ub', 'guess', 'par', 'stats', 'status', 'obj', ...
                  'horizon', 'verbosity', 'isQP', 'name', 'timelimit', ...
                  'maxiter', 'solveroptions'};
    end%function
    
    function [conlb_, conub_] = getconbounds(self)
        % Returns bounds for constraints, accounting for any truncated horizon.
        conlb_ = self.conlb;
        conub_ = self.conub;
        if self.horizon ~= self.maxhorizon
            n = self.maxhorizon - self.horizon;
            fields = {'measurements', 'model'};
            for i = 1:length(fields)
                f = fields{i};
                if isfield(conlb_, f)
                    conlb_.(f) = mpctools.setfinalvalues(conlb_.(f), n, -inf());
                end
                if isfield(conub_, f)
                    conub_.(f) = mpctools.setfinalvalues(conub_.(f), n, inf());
                end
            end
        end
    end%function
    
    function [lb_, ub_] = getvarbounds(self)
        % Returns bounds for variables.
        %
        % Intended to be overridden in subclasses to allow for special handling
        % of truncated horizons.
        lb_ = self.lb;
        ub_ = self.ub;
    end%function
end%methods

end%classdef
