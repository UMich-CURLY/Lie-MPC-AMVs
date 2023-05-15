classdef MHESolver < mpctools.ControlSolver
% Class for holding an NLP solver for MHE problems. Note that the
% user should not create instances of this class directly (i.e., using the
% constructor), but should instead use the `nmhe` function.
%
% In addition to all the properties from the `ControlSolver` class, this object
% also contains the following extra fields:
%
% * `history`
%
% > Array struct that holds a history of estimates for $x$, $w$, and $v$, as
% well as values of $y$ and $u$. `self.history(t)` gives the estimate from
% `t` time periods ago. The current solution is stored to the history whenever
% `self.saveestimate()` is called.
%
% > Keyword argument 'Nhistory' controls the number of past estimates to
% include in `history`. The default value is one more than the horizon of the
% MHE problem.
%
% * `priorupdate` : String
%
% > String indicating which type of prior update (if any) is being used.
% See documentation of `nmhe()` for more details.
%
% * `fix_truncated_x` : Logical [`false`]
%
% > Whether to fix the $x$ variables that are outside of the current horizon.
% Usually this is unnecessary, but it can help with
% `"Restoration_Failed"` issues when the horizon is truncated. Note that you may
% need to provide a feasible guess for $x$ if `fix_truncated_x` is true.
properties (Access=public)
    fix_truncated_x = false();
end%properties

properties (GetAccess=public, SetAccess=protected)
    history = struct();
    priorupdate = 'none';
end%properties

properties (GetAccess=protected, SetAccess=protected)
    priorfuncs = struct();
end%properties

methods (Access=public)
    function self = MHESolver(varargin)
        % self = MHESolver(...)
        %
        % Creates the solver object.
        persistent parser
        if isempty(parser)
            parser = mpctools.ArgumentParser();
            parser.add('Nhistory', [], {'positive', 'scalar', 'int'});
            parser.add('priorupdate', 'none', 'str');
            parser.add('mhemodel', [], 'struct');
            parser.add('fix_truncated_x', false(), 'bool');
            parser.add('**kwargs');
        end
        [args, superargs] = parser.parse(varargin{:});
        
        % Call superclass constructor.
        self = self@mpctools.ControlSolver('**', superargs);
        
        % Preallocate the history struct.
        self.resethistory(args.Nhistory);
        
        % Calculate information for prior update.
        availableupdates = {'none', 'filtering', 'smoothing', 'hybrid'};
        if ~ismember(args.priorupdate, availableupdates)
            error('Unknown prior update method ''%s''!', args.priorupdate);
        end
        self.priorupdate = args.priorupdate;
        self.initpriorfuncs(args.mhemodel);
        self.fix_truncated_x = args.fix_truncated_x;
    end%function
    
    function newmeasurement(self, y, u, x0bar)
        % `self.newmeasurement(y, [u], [x0bar])`
        %
        % Shifts `par.y` to include the new measurement `y` at the current time.
        % Also shifts `par.u` if `u` is given. If `x0bar` is given, its value is
        % simply overwritten.
        %
        % In contrast to `cyclepar()`, this function accounts for a shortened
        % horizon, and it only removes old data if the horizon is full.
        narginchk(2, 4);
        
        if ~isfield(self.par, 'y')
            error('y is missing from par!');
        end
        if self.horizon < self.maxhorizon
            self.par.y(:,self.horizon + 1) = y(:);
        else
            self.par.y = [self.par.y(:,2:end), y(:)];
        end

        if nargin() > 2 && ~isempty(u)
            if ~isfield(self.par, 'u')
                error('u is missing from par!');
            end
            if self.horizon < self.maxhorizon
                self.par.u(:,self.horizon + 1) = u(:);
            else
                self.par.u = [self.par.u(:,2:end), u(:)];
            end
        end

        if nargin() > 3 && ~isempty(x0bar)
            self.par.x0bar = x0bar;
        end
    end%function
    
    function saveestimate(self, updateguess)
        % `self.saveestimate([updateguess=true])`
        %
        % Saves the current estimates to the history struct. Also updates the
        % prior weight if `self.priorupdate` has been specified.
        %
        % If `updateguess` is `true` and the horizon is full, then the solver's
        % guess is also updated.
        narginchk(1, 2);
        if nargin() < 2
            updateguess = true();
        end
        
        % Get current value of history.
        thish = struct();
        fields = fieldnames(self.history);
        for i = 1:length(fields)
            f = fields{i};
            switch f
            case {'xhat', 'what', 'vhat'}
                thish.(f) = self.var.(f(1));
            case 'yhat'
                thish.(f) = self.par.y - self.var.v;
            case {'y', 'u'}
                thish.(f) = self.par.(f);
            otherwise
                warning('Unknown field in history: ''%s''', f);
                thish.(f) = [];
            end
        end
        
        % Truncate variables if necessary.
        Ntrunc = self.maxhorizon - self.horizon;
        if Ntrunc > 0
            thish = structfun(@(x) x(:,1:(end - Ntrunc)), thish, ...
                              'UniformOutput', false());
        end
        
        % Add new history to front.
        self.history = [thish; self.history(1:(end - 1))];
        
        % Update prior if horizon is full.
        if Ntrunc == 0
            switch self.priorupdate
            case 'filtering'
                self.filteringupdate();
            case 'smoothing'
                self.smoothingupdate();
            case 'hybrid'
                self.hybridupdate();
            end
        end
        
        % Update guess.
        if updateguess
            self.saveguess(self.var, double(Ntrunc == 0));
        end
    end%function
    
    function reset(self)
        % `self.reset()`
        %
        % Resets the object to its initial state.
        narginchk(1, 1);
        self.resethistory(length(self.history));
        reset@mpctools.ControlSolver(self);
    end%function
end%methods

methods (Access=protected)
    function resethistory(self, N)
        % self.resethistory([Nhistory])
        %
        % Resets the `history` and `t` fields of the object. Optional argument
        % `Nhistory` decides how many time points to save in the history
        % object.
        narginchk(1, 2)
        if nargin() < 2 || isempty(N)
            N = self.maxhorizon + 1;
        end
        self.history = struct();
        h = cell(N, 1);
        historyfields = {'xhat', h, 'what', h, 'vhat', h, 'yhat', h, 'y', h};
        if isfield(self.par, 'u')
            historyfields = [historyfields, {'u', h}];
        end
        self.history = struct(historyfields{:});
    end%function
    
    function initpriorfuncs(self, model)
        % Initializes the functions used for the prior update.
        self.priorfuncs = struct();
        if isequal(self.priorupdate, 'none')
            return
        elseif isempty(model)
            error('Must provide model functions to use prior update!');
        end
        expr = struct('x', self.nlp.x, 'p', self.nlp.p);
        
        % Get one-step prediction function and linearizations.
        sym = mpctools.structupdate(self.parsym, self.varsym);
        shootingargs = {model.f.func, model.f.args, model.h.func, model.h.args};
        onestep = mpctools.singleshooting(sym, shootingargs{:}, 1:2);
        
        measure = expr;
        measure.h = onestep.y{1};
        measure.C = measure.h.jacobian(sym.x{1});
        measure.H = measure.h.jacobian(sym.v{1});
        self.priorfuncs.measure = ...
            casadi.Function('measure', measure, {'x', 'p'}, {'h', 'C', 'H'});
        
        predict = expr;
        predict.f = onestep.x{2};
        predict.A = predict.f.jacobian(sym.x{1});
        predict.G = predict.f.jacobian(sym.w{1});
        self.priorfuncs.predict = ...
            casadi.Function('predict', predict, {'x', 'p'}, {'f', 'A', 'G'});
        
        % Get linearized objective function.
        obj = expr;
        obj.l = self.objsym;
        obj.Pinv = 0.5*obj.l.hessian(sym.x{1});
        obj.Qinv = 0.5*obj.l.hessian(sym.w{1});
        obj.Rinv = 0.5*obj.l.hessian(sym.v{1});
        self.priorfuncs.obj = ...
            casadi.Function('obj', obj, {'x', 'p'}, ...
                            {'l', 'Pinv', 'Qinv', 'Rinv'});
        
        % Compute some extra functions.
        switch self.priorupdate
        case {'filtering', 'hybrid'}
            % Nothing left to do.
        case 'smoothing'
            shooting = mpctools.singleshooting(sym, shootingargs{:}, ...
                                               2:length(sym.x));
            
            W = cell(2, length(shooting.v) - 1);
            W(1,:) = shooting.v(2:end);
            W(2,1:(end - 1)) = shooting.w(2:end);
            W = vertcat(W{1:(end-1)}); % Interleaved w and v.
            
            expr.Y = vertcat(shooting.y{2:end});
            expr.O = expr.Y.jacobian(shooting.x{2});
            expr.G = expr.Y.jacobian(W);
            
            expr.Qinv = 0.5*self.objsym.hessian(W);
            
            self.priorfuncs.shooting = ...
                casadi.Function('shooting', expr, {'x', 'p'}, ...
                                {'Y', 'O', 'G', 'Qinv'});
        end
    end%function
    
    function [xhatm, Pm] = ekfstep(self, xhatm)
        % [xhatm, Pm] = self.ekfstep(xhatm)
        %
        % Returns updated xhatm and Pm based on one step of the EKF applied to
        % the initial state estimate. For time-varying f and h, uses t = 0 for
        % f and t = 1 for h.
        %
        % The function takes xhat(t | t) and simulates without disturbances to
        % compute xhat(t + 1 | t). The covariance for xhat(t + 1 | t) is also
        % calculated.
        narginchk(2, 2);
        p = self.map.p.vec(self.par);
        
        thisvar = self.var;
        thisvar.x(:,1) = xhatm;
        thisvar.v(:,1) = 0;
        thisvar.w(:,1) = 0;
        
        mats = self.priorfuncs.predict('x', self.map.x.vec(thisvar), 'p', p);
        xhatm = full(mats.f);
        A = full(mats.A);
        G = full(mats.G);
        
        mats = self.priorfuncs.measure('x', self.map.x.vec(thisvar), 'p', p);
        C = full(mats.C);
        H = full(mats.H);
        
        mats = self.priorfuncs.obj('x', self.map.x.vec(thisvar), 'p', p);
        P = mpctools.spdinv(full(mats.Pinv));
        Q = mpctools.spdinv(full(mats.Qinv));
        R = mpctools.spdinv(full(mats.Rinv));
        
        Pm = G*Q*G' + A*P*A' - (A*P*C')/(C*P*C' + H*R*H')*(C*P*A');
    end%function
    
    function filteringupdate(self)
        % Update prior covariance using filtering update.
        narginchk(1, 1);
        
        if length(self.history) <= self.horizon
            error('Not enough history is present to use the filtering update!');
        end
        xhat = self.history(self.horizon + 1).xhat(:,end);
        [self.par.x0bar, Pm] = self.ekfstep(xhat);
        self.par.Pinv = mpctools.spdinv(Pm);
    end%function
    
    function smoothingupdate(self)
        % Update prior covariance using smoothing update.
        narginchk(1, 1);
        
        xhat = self.history(1).xhat(:,1);
        [xhatm, Pm] = self.ekfstep(xhat);
        
        thisvar = self.var;
        thisvar.w(:) = 0;
        thisvar.v(:) = 0;
        thisvar.x(:,2) = xhatm;
        
        mats = self.priorfuncs.shooting('x', self.map.x.vec(thisvar), ...
                                        'p', self.map.p.vec(self.par));
        
        G = full(mats.G);
        O = full(mats.O);
        Q = mpctools.spdinv(full(mats.Qinv));
        Y = reshape(self.par.y(:,2:end), [], 1);
        Ybar = full(mats.Y);
        E = Y - Ybar;
        
        L = (Pm*O')/(G*Q*G' + O*Pm*O');
        
        mxy = self.history(1).xhat(:,2);
        Pxy = Pm - L*O*Pm;
        Pxyinv = mpctools.spdinv(Pxy);
        
        Pyx = G*Q*G';
        Pyxinv = mpctools.spdinv(Pyx);
        
        % Sum the quadratic and linear terms, and then factor.
        Hess = Pxyinv - O'*Pyxinv*O;
        grad = Pxyinv*mxy - O'*Pyxinv*(E + O*xhatm);
        
        self.par.Pinv = Hess;
        self.par.x0bar = Hess\grad;
    end%function
    
    function hybridupdate(self)
        % Update Pinv using filtering covariance but x0bar using smoothing
        % estimate.
        narginchk(1, 1);
        
        xhat = self.history(1).xhat(:,1);
        [self.par.x0bar, Pm] = self.ekfstep(xhat);
        self.par.Pinv = mpctools.spdinv(Pm);
    end%function
    
    function [lb_, ub_] = getvarbounds(self)
        % Returns bounds for variables, fixing any excluded x, y, and z
        % variables.
        lb_ = self.lb;
        ub_ = self.ub;
        Ntrunc = self.maxhorizon - self.horizon;
        if Ntrunc > 0
            fixedvars = {'x', 'y', 'z'};
            fixedvars = [fixedvars, cellfun(@(v) [v, 'c'], fixedvars, ...
                                            'UniformOutput', false())];
            for i = 1:length(fixedvars)
                v = fixedvars{i};
                if isfield(self.varsym, v)
                    lb_.(v) = mpctools.setfinalvalues(lb_.(v), Ntrunc, ...
                                                      self.guess.(v));
                    ub_.(v) = mpctools.setfinalvalues(ub_.(v), Ntrunc, ...
                                                      self.guess.(v));
                end
            end
        end
    end%function
end%methods

end%classdef

