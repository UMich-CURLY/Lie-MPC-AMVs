classdef KalmanFilter < handle
% Class for (offset-free) Kalman Filter.
%
% Public attributes are as follows:
%
% * `N`
%
% > struct of system sizes. Fields are `x`, `u`, `y`, `d`, `w`, and `v`.
%
% * `A`
% * `B`
% * `C`
% * `D`
% * `Ad`
% * `Bd`
% * `Cd`
% * `Gx`
% * `Gd`
% * `Qw`
% * `Rv`
%
% > Matrices that define the system and disturbance model.
%
%
% * `xhat`
% * `dhat`
%
% > Current estimates of state and integrating disturbances after the current
% measurement has been considered. These are updated via calls to `filter()`.
%
% * `xhatm`
% * `dhatm`
%
% > Current estimates of state and integrating disturbances before the current
% measurement has been considered. These are updated via calls to `predict()`.
%
% * `Lx`
% * `Ld`
%
% > Kalman Filter gains for the states and integrating disturbances.
%
% * `Txd`
% * `Txy`
% * `Tud`
% * `Tuy`
%
% > Matrices for the steady-state target calculation. Multiply either d or ysp
% and return either xtarg or utarg.
properties (GetAccess=public, SetAccess=public)
    xhatm;
    dhatm;
    xhat;
    dhat;
    u;
end%properties

properties (GetAccess=public, SetAccess=private)
    N;
    A;
    B;
    C;
    D;
    Ad;
    Bd;
    Cd;
    Qw;
    Rv;
    Gx;
    Gd;
    Lx;
    Ld;
    H;
    Txd;
    Txr;
    Tud;
    Tur;
    xs;
    us;
    ys;
end%properties

properties (GetAccess=private, SetAccess=private)
    sizes = struct('A', 'xx', 'B', 'xu', 'C', 'yx', 'D', 'yu', ...
                   'Bd', 'xd', 'Cd', 'yd', 'Qw', 'ww', 'Rv', 'vv', ...
                   'Gx', 'xw', 'Gd', 'dw', 'Ad', 'dd', 'H', 'ry', ...
                   'xs', 'x', 'ys', 'y', 'us', 'u');
end%properties

methods (Access=public)
    function self = KalmanFilter(varargin)
        % `self = KalmanFilter(...)`
        %
        % Creates a KalmanFilter instance.
        %
        % Arguments are as follows:
        %
        % * `A` : Matrix
        % * `B` : Matrix
        % * `C` : Matrix
        %
        % > Matrices to define nominal model.
        %
        % * `distmodel` : string
        %
        % > Either 'input', 'output', 'custom', or 'none' to specify what kind
        % of disturbance model to use.
        %
        % > Default is 'custom'.
        %
        % * `Bd` : Matrix
        % * `Cd' : Matrix
        %
        % > Matrices to use if `distmodel` is set to 'custom'.
        %
        % * `Qw` : Matrix
        % * `Rv` : Matrix
        %
        % > Covariance matrices to use for the Kalman Filter. Both default to
        % identity matrix of appropriate size.
        %
        % * `Gx` : Matrix
        % * `Gd` : Matrix
        %
        % > State noise shaping matrices. System evolution is assumed to be
        %
        %     x^+ = A*x + B*u + Bd*d + Gx*w
        %     d^+ = Ad*d + Gd*w
        %     y = C*x + Cd*d + D*u + v
        %
        % > If not specified, Gx and Gd are chosen as appropriate sections of
        % the identity matrix.
        %
        % * `D` : Matrix
        %
        % > Direct input-to-output feedthrough matrix. Note that you must change
        % the `u` attribute of the object to the initial value of u.
        %
        % * `Ad` : Matrix
        %
        % > Matrix to use for integrating disturbance evolution. Default is
        % identity.
        %
        % * `ifundetectable` : string
        %
        % > Specifies what to do if the augmented disturbance model is
        % undetectable. Choices are 'error', 'warning', or 'ignore'. Default is
        % 'warning'.
        %
        % * `contvars` : Row Vector
        %
        % > Gives the indices of outputs that are controlled. Only used in the
        % steady-state target problem. Defaults to 1:Nu.
        %
        % * `H` : Matrix
        %
        % > Calculates `r = H*y` with `r` the given setpoint. Note that r should
        % be smaller than `u`, or else the steady-state target problem may be
        % infeasible.
        %
        % > Note that specifying this option overrides the choice of `contvars`.
        %
        % * `xs`, `us`, `us` : Column Vectors
        %
        % > Steady-state values if using positional variables. All default to
        % zero.
        persistent parser
        if isempty(parser)
            parser = mpctools.ArgumentParser();
            parser.add('A', 'required', 'matrix');
            parser.add('B', 'required', 'matrix');
            parser.add('C', 'required', 'matrix');
            parser.add('distmodel', 'custom', 'str');
            parser.add('Bd', [], 'matrix');
            parser.add('Cd', [], 'matrix');
            parser.add('Qw', [], 'matrix');
            parser.add('Rv', [], 'matrix');
            parser.add('Gx', [], 'matrix');
            parser.add('Gd', [], 'matrix');
            parser.add('D', [], 'matrix');
            parser.add('Ad', [], 'matrix');
            parser.add('ifundetectable', 'warning', 'str');
            parser.add('contvars', [], {'row', 'int', 'positive'});
            parser.add('H', [], 'matrix');
            parser.add('xs', [], 'column');
            parser.add('ys', [], 'column');
            parser.add('us', [], 'column');
        end
        [args, ~, isdefault] = parser.parse(varargin{:});
        
        % Get sizes from A, B, and C.
        self.N = struct('x', size(args.A, 1), 'u', size(args.B, 2), ...
                        'y', size(args.C, 1));
        self.N.v = self.N.y; % This guy is fixed.
        
        % Choose disturbance model.
        switch args.distmodel
        case 'output'
            self.N.d = self.N.y;
            isdefault.Bd = true();
            isdefault.Cd = false();
            args.Cd = eye(self.N.d);
        case 'input'
            self.N.d = self.N.u;
            isdefault.Bd = false();
            isdefault.Cd = true();
            args.Bd = args.B;
        case 'none'
            isdefault.Bd = true();
            isdefault.Cd = true();
        case 'custom'
            if isdefault.Bd && isdefault.Cd
                error('At least one of Bd and Cd must be given!');
            end
            if isdefault.Bd
                self.N.d = size(args.Cd, 2);
            else
                self.N.d = size(args.Bd, 2);
            end
        end
        
        % Check sizes and fill in defaults.
        args = self.checksizes(args, isdefault);
        self.A = args.A;
        self.B = args.B;
        self.C = args.C;
        self.D = args.D;
        self.Ad = args.Ad;
        self.Bd = args.Bd;
        self.Cd = args.Cd;
        self.Qw = args.Qw;
        self.Rv = args.Rv;
        self.Gx = args.Gx;
        self.Gd = args.Gd;
        self.H = args.H;
        self.xs = args.xs;
        self.ys = args.ys;
        self.us = args.us;
        
        % Compute Kalman filter.
        Aaug = [self.A, self.Bd; zeros(self.N.d, self.N.x), self.Ad];
        Baug = [self.B; zeros(self.N.d, self.N.u)];
        Caug = [self.C, self.Cd];
        [isdetec, nbad] = mpctools.isdetectable(Aaug, Baug, Caug);
        if ~isdetec
            msg = sprintf('Augmented system has %d undetectable modes!', nbad);
        elseif self.N.d < self.N.y
            msg = 'Not enough disturbances added. System may have offset!';
        else
            msg = '';
        end
        if ~isempty(msg)
            switch args.ifundetectable
            case 'ignore'
                % Pass.
            case {'warning', 'warn'}
                warning(msg);
            case 'error'
                error(msg)
            otherwise
                error('Invalid choice for ifundetectable!');
            end
        end
        Gaug = [self.Gx; self.Gd];
        L = dlqe(Aaug, Gaug, Caug, self.Qw, self.Rv);
        self.Lx = L(1:self.N.x,:);
        self.Ld = L((self.N.x + 1):end,:);
        
        % Compute steady-state target matrices.
        Ttarg = [eye(self.N.x) - self.A, -self.B; self.H*[self.C, self.D]];
        Ttarg = Ttarg\[self.Bd, zeros(self.N.x, self.N.r); ...
                       -self.H*self.Cd, eye(self.N.r)];
        self.Txd = Ttarg(1:self.N.x,1:self.N.d);
        self.Txr = Ttarg(1:self.N.x,(self.N.d + 1):end);
        self.Tud = Ttarg((self.N.x + 1):end,1:self.N.d);
        self.Tur = Ttarg((self.N.x + 1):end,(self.N.d + 1):end);
        
        % Initial guesses for xhat and dhat.
        self.xhat = zeros(self.N.x, 1);
        self.xhatm = self.xhat;
        self.dhat = zeros(self.N.d, 1);
        self.dhatm = self.dhat;
        self.u = zeros(self.N.u, 1);
    end%function
    
    function [xhat_, dhat_] = filter(self, y, xhatm_, dhatm_, u_)
        % `[xhat, dhat] = self.filter(y, [xhatm], [dhatm], [u])`
        %
        % Performs the Kalman Filter update using the current measurement y.
        % New values are stored to `self.xhat` and `self.dhat`.
        %
        % If not given, `xhatm` and `dhatm` are taken from the current values of
        % the object. `u` is handled similarly, but it only matters if `self.D`
        % is nonzero.
        narginchk(2, 5);
        if nargin() < 3
            xhatm_ = self.xhatm;
        end
        if nargin() < 4
            dhatm_ = self.dhatm;
        end
        if nargin() < 5
            u_ = self.u;
        end
        
        e = y - self.ys - (self.C*(xhatm_ - self.xs) + self.Cd*dhatm_ ...
                           + self.D*(u_ - self.us));
        xhat_ = xhatm_ + self.Lx*e;
        dhat_ = dhatm_ + self.Ld*e;
        
        self.xhat = xhat_;
        self.dhat = dhat_;
    end%function
    
    function [xhatm_, dhatm_] = predict(self, u_, xhat_, dhat_)
        % `[xhatm, dhatm] = self.predict(u, [xhat], [dhat])`
        %
        % Performs the Kalman Filter prediction using the current control
        % action `u`. New values are stored to `self.xhatm` and `self.dhatm`.
        %
        % If not given, `xhat` and `dhat` are taken from the current values of
        % the object.
        narginchk(2, 4);
        if nargin() < 3
            xhat_ = self.xhat;
        end
        if nargin() < 4
            dhat_ = self.dhat;
        end
        
        xhatm_ = self.A*(xhat_ - self.xs) + self.B*(u_ - self.us) ...
                 + self.Bd*dhat_ + self.xs;
        dhatm_ = self.Ad*dhat_;
        
        self.xhatm = xhatm_;
        self.dhatm = dhatm_;
        self.u = u_;
    end%function
    
    function [xtarg, utarg] = target(self, varargin)
        % `[xtarg, utarg] = self.target(ysp, [dhat], [rsp])`
        %
        % Calculates the steady-state target values of x and u using the given
        % value of `ysp` or `rsp`.
        %
        % If `rsp` is given, it is used directly. Otherwise, it is calculated as
        % `rsp = self.H*(ysp - self.ys)`.
        persistent parser
        if isempty(parser)
            parser = mpctools.ArgumentParser();
            parser.add('ysp', [], 'col');
            parser.add('dhat', [], 'col');
            parser.add('rsp', [], 'col');
        end
        [args, ~, isdefault] = parser.parse(varargin{:});
        if isdefault.rsp
            if isdefault.ysp
                error('Either ysp or rsp must be given!');
            end
            args.rsp = self.H*(args.ysp - self.ys);
        end
        if isdefault.dhat
            args.dhat = self.dhat;
        end
        xtarg = self.xs + self.Txd*args.dhat + self.Txr*args.rsp;
        utarg = self.us + self.Tud*args.dhat + self.Tur*args.rsp;
    end%function
end%methods

methods (Access=private)
    function args = checksizes(self, args, isdefault)
        % Checks the sizes of given arguments the class constructor,
        % saving sizes and filling in appropriate defaults.
        narginchk(3, 3);
        names = fieldnames(self.sizes);
        
        % Loop through arguments to get sizes.
        for i = 1:length(names)
            n = names{i};
            if ~isdefault.(n)
                sizesyms = self.sizes.(n);
                for j = 1:length(sizesyms)
                    s = sizesyms(j);
                    if ~isfield(self.N, s)
                        self.N.(s) = size(args.(n), j);
                    end
                end
            end
        end
        if ~isfield(self.N, 'd')
            self.N.d = 0;
        end
        if ~isfield(self.N, 'w')
            self.N.w = self.N.x + self.N.d;
        end
        
        % Set H as special case.
        if isdefault.H
            isdefault.H = false();
            if isdefault.contvars
                args.contvars = 1:self.N.u;
            end
            args.H = eye(self.N.y);
            args.H = args.H(args.contvars,:);
        end
        self.N.r = size(args.H, 1);
        
        % Set other default matrices.
        for i = 1:length(names)
            n = names{i};
            if isdefault.(n)
                args.(n) = self.getdefault(n, args.(n));
            end
        end
        
        % Check sizes of all matrices.
        okayrows = true(length(names), 1);
        badstrs = cell(length(names), 1);
        for i = 1:length(names)
            n = names{i};
            issize = size(args.(n));
            shouldbesize = self.getsize(n);
            if ~isequal(issize, shouldbesize)
                okayrows(i) = false();
                badstrs{i} = sprintf('    %s is %s but should be %s', n, ...
                                     mpctools.row2str(issize), ...
                                     mpctools.row2str(shouldbesize));
            end
        end
        if ~all(okayrows)
            badstrs = badstrs(~okayrows);
            error('Incorrect sizes:\n%s', sprintf('%s\n', badstrs{:}));
        end
    end%function
    
    function M = getdefault(self, name, scale)
        % Returns default matrix for argument. Certain matrices are scaled by
        % scale if it is a scalar.
        if ~isscalar(scale)
            scale = 1;
        end
        narginchk(3, 3);
        switch name
        case 'Gx'
            M = eye(self.N.x, self.N.w);
        case 'Gd'
            M = [zeros(self.N.d, self.N.w - self.N.d), eye(self.N.d)];
        case 'Qw'
            M = scale*eye(self.N.w);
        case 'Rv'
            M = scale*eye(self.N.v);
        case 'Ad'
            M = eye(self.N.d);
        otherwise
            % Default case is appropriately-sized zeros.
            M = zeros(self.getsize(name));
        end
    end%function
    
    function sz = getsize(self, name)
        % Returns the correct numeric size for the given matrix.
        sizesym = self.sizes.(name);
        sz = zeros(size(sizesym));
        for i = 1:length(sizesym)
            sz(i) = self.N.(sizesym(i));
        end
        if length(sz) < 2
            sz = [sz, 1];
        end
    end%function
end%methods

end%classdef
