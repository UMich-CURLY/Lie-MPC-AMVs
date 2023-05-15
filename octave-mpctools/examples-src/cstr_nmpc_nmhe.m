% Control of a nonlinear CSTR using nonlinear MPC and MHE.
pkg('load', 'control');
mpc = import_mpctools();

% Choose disturbance model and whether to use nonlinear MPC.
disturbancemodel = 'Good'; % Choose 'No', 'Good', 'Offset', 'Undetectable', or 'Exact'
nonlinearMPC = true();

% Sizes for the nonlinear system.
Delta = 1;
Nx = 3;
Nu = 2;
Ny = Nx;
Np = 1;
Nw = Nx;
Nv = Ny;
small = 1e-5; % Small number.

% Parameters.
pars = struct();
pars.T0 = 350; % K
pars.c0 =  1;  % kmol/m^3
pars.r = 0.219; % m
pars.k0 = 7.2e10; % min^-1
pars.E =  8750; % K
pars.U =  54.94; % kJ/(min m^2 K)
pars.rho = 1e3;   % kg/m^3
pars.Cp = 0.239;  % kJ/(kg K)
pars.DeltaH = -5e4; % kJ/kmol
pars.A = pi()*pars.r.^2;
pars.rhoCp = pars.rho*pars.Cp;

function rhs = cstrode(x, u, p, pars)
    % Nonlinear ODE model for reactor.
    c = x(1);
    T = x(2);
    h = x(3) + eps(); % Avoids division by zero.
    
    Tc = u(1);
    F = u(2);
    
    F0 = p(1);
    
    k = pars.k0*exp(-pars.E/T);
    rate = k*c;
    
    dcdt = F0*(pars.c0 - c)/(pars.A*h) - rate;
    dTdt = F0*(pars.T0 - T)/(pars.A*h) ...
           - pars.DeltaH/pars.rhoCp*rate ...
           + 2*pars.U/(pars.r*pars.rhoCp)*(Tc - T); 
    dhdt = (F0 - F)/pars.A;
    
    rhs = [dcdt; dTdt; dhdt];
end%function
ode = @(x, u, p) cstrode(x, u, p, pars);

cstrsim = mpc.getCasadiIntegrator(ode, Delta, [Nx, Nu, Np], ...
                                  {'x', 'u', 'p'}, {'cstr'});

% Steady-state values.
cs = 0.878; % kmol/m^3
Ts = 324.5; % K
hs = 0.659; % m
Fs = 0.1; % m^3/min
Tcs = 300; % K
F0s = 0.1; % m^3/min

xs = [cs; Ts; hs];
us = [Tcs; Fs];
ps = [F0s];

CVs = [1, 3]; % Control concentration and height.

xlb = [0.01; 300; 0.01];
xub = [10; 350; 10];

% Simulate a few steps so we actually get dx/dt = 0 at steady state.
for i = 1:10
    xs = full(cstrsim(xs, us, ps));
end
cs = xs(1);
Ts = xs(2);
hs = xs(3);
C = eye(Ny, Nx);
ys = C*xs;

% Bounds.
umax = [0.05*Tcs; 0.5*Fs];
ulb = us - umax;
uub = us + umax;

% Linearize the model at this point.
[A, B, Bp] = mpc.getLinearizedModel(ode, {xs, us, ps}, 'Delta', Delta, ...
                                    'deal', true());
flin = @(x, u, p) A*(x - xs) + B*(u - us) + Bp*(p - ps) + xs;

% Choose disturbance model.
fprintf('Choosing %s disturbance model.\n', disturbancemodel);
switch disturbancemodel
case 'Good'
    Nd = 3;
    Bd = [0, 0, 0; 0, 0, 1];
    Cd = [1 0 0; 0 0 0; 0 1 0];
    Bpd = zeros(Np, Nd);
    Qwd = diag([small, small, 1]);
case 'Offset'
    Nd = 2;
    Bd = zeros(Nu, Nd);
    Cd = [1 0; 0 0; 0 1];
    Bpd = zeros(Np, Nd);
    Qwd = diag([small, 1]);
case 'Undetectable'
    Nd = 3;
    Bd = zeros(Nu, Nd);
    Cd = eye(Ny, Nd);
    Bpd = zeros(Np, Nd);
    Qwd = diag([small, small, 1]);
case 'No'
    Nd = 0;
    Bd = zeros(Nu, Nd);
    Cd = zeros(Ny, Nd);
    Bpd = zeros(Np, Nd);
    Qwd = eye(Nd);
case 'Exact'
    Nd = 1;
    Bd = zeros(Nu, Nd);
    Cd = zeros(Ny, Nd);
    Bpd = eye(Np, Nd);
    Qwd = 1;
otherwise
    error('Unknown choice for disturbance model: %s', disturbancemodel);
end

% Get Kalman Filter weight as a prior.
Qw = small*eye(Nx);
if isequal(disturbancemodel, 'No')
    Qw(end,end) = 1;
end
Rv = small*diag(ys.^2);

Aaug = [A, B*Bd + Bp*Bpd; zeros(Nd, Nx), eye(Nd)];
Baug = [B; zeros(Nd, Nu)];
Gaug = eye(Nx + Nd);
Caug = [eye(Ny, Nx), Cd];
Qaug = blkdiag(Qw, Qwd);
Raug = Rv;

detec = rank([eye(Nx + Nd) - Aaug; Caug]);
if detec < Nx + Nd
    fprintf(' * Augmented system is not detectable!\n')
    Aaug = Aaug - 1e-6*eye(size(Aaug)); % Help out dlqe a bit.
end
[Lkf, ~, Pkf] = dlqe(Aaug, Gaug, Caug, Qaug, Raug);

% Make functions for controller, target finder, and estimator.
fprintf('Building controllers.\n');

if nonlinearMPC
    model = ode;
    discretization = {'rk4', true(), 'Delta', Delta, 'M', 2};
else
    model = flin;
    discretization = {};
end
modelargs = {@(x, u, Bd, Bpd, d) model(x, u + Bd*d, ps + Bpd*d), ...
             {Nx, Nu, [Nu, Nd], [Np, Nd], Nd}, {'x', 'u', 'Bd', 'Bpd', 'd'}};
f = mpc.getCasadiFunc(modelargs{:}, discretization{:}, ...
                      'funcname', 'fcontroller');
fsstarg = mpc.getCasadiFunc(modelargs{:}, discretization{:}, ...
                            'funcname', 'fsstarg');

% Make stage costs.
function cost = stagecost(x, u, xsp, usp, Q, R)
    dx = x - xsp;
    du = u - usp;
    cost = dx'*Q*dx + du'*R*du;
end%function
lcontroller = mpc.getCasadiFunc(@stagecost, ...
                                {Nx, Nu, Nx, Nu, [Nx, Nx], [Nu, Nu]}, ...
                                {'x', 'u', 'xsp', 'usp', 'Q', 'R'}, ...
                                'funcname', 'l');

Vf = mpc.getCasadiFunc(@(x, xsp, P) stagecost(x, 0, xsp, 0, P, 0), ...
                       {Nx, Nx, [Nx, Nx]}, {'x', 'xsp', 'P'}, 'funcname', 'Vf');

lmhe = mpc.getCasadiFunc(@(w, v, Dd, Qinv, Rinv, Qdinv) ...
                           w'*Qinv*w + v'*Rinv*v + Dd'*Qdinv*Dd, ...
                         {Nw, Nv, Nd, [Nw, Nw], [Nv, Nv], [Nd, Nd]}, ...
                         {'w', 'v', 'Dd', 'Qinv', 'Rinv', 'Qdinv'}, ...
                         'funcname', 'l');

function cost = priorcost(x, d, xhat, dhat, Pinv)
    z = [x; d];
    zhat = [xhat; dhat];
    dz = z - zhat;
    cost = dz'*Pinv*dz;
end%function
lprior = mpc.getCasadiFunc(@priorcost, ...
                           {Nx, Nd, Nx, Nd, [Nx + Nd, Nx + Nd]}, ...
                           {'x', 'd', 'xbar', 'dbar', 'Pinv'}, ...
                           'funcname', 'lx');

h = mpc.getCasadiFunc(@(x, Cd, d) x + Cd*d, {Nx, [Ny, Nd], Nd}, ...
                      {'x', 'Cd', 'd'}, 'funcname', 'h');

% Assemble controller, target finder, and estimator.
Ncontroller = 5;
Nmhe = 5;
Q = diag(xs.^-2);
R = diag(us.^-2);
[~, P] = dlqr(A, B, Q, R);

N = struct('x', Nx, 'u', Nu, 't', Ncontroller);
guess = struct('x', repmat(xs, 1, N.t + 1), 'u', repmat(us, 1, N.t));
lb = struct('x', xlb, 'u', ulb);
ub = struct('x', xub, 'u', uub);
par = struct('xsp', xs, 'usp', us, 'Q', Q, 'R', R, 'P', P, ...
             'Bd', Bd, 'Bpd', Bpd, 'd', zeros(Nd, 1));
controller = mpc.nmpc('f', f, 'l', lcontroller, 'Vf', Vf, 'N', N, ...
                      'lb', lb, 'ub', ub, 'guess', guess, 'par', par);

N = struct('x', Nx, 'u', Nu, 'y', Ny);
guess = struct('x', xs, 'u', us, 'y', ys);
lb = struct('x', xlb, 'u', ulb);
ub = struct('x', xub, 'u', uub);
par = struct('Bd', Bd, 'Bpd', Bpd, 'Cd', Cd, 'd', zeros(Nd, 1));
sstarg = mpc.sstarg('f', fsstarg, 'h', h, 'N', N, ...
                    'lb', lb, 'ub', ub, 'guess', guess, 'par', par);

N = struct('x', Nx, 'd', Nd, 'u', Nu, 'y', Ny, 't', Nmhe);
guess = struct('x', repmat(xs, 1, N.t + 1));
lb = struct('x', xlb);
ub = struct('x', xub);
par = struct('y', repmat(ys, 1, N.t + 1), 'u', repmat(us, 1, N.t), ...
             'xbar', xs, 'dbar', zeros(Nd, 1), ...
             'Qinv', mpc.spdinv(Qw), 'Rinv', mpc.spdinv(Rv), ...
             'Qdinv', mpc.spdinv(Qwd), 'Pinv', mpc.spdinv(Pkf), ...
             'Bd', Bd, 'Cd', Cd, 'Bpd', Bpd);
mhe = mpc.nmhe('f', f, 'h', h, 'l', lmhe, 'lx', lprior, 'N', N, ...
               'guess', guess, 'lb', lb, 'ub', ub, 'par', par, ...
               'wadditive', true());

% Closed-loop simulation.
Nsim = 40;
x = NaN(Nx, Nsim + 1);
x(:,1) = xs;
y = NaN(Ny, Nsim + 1);
u = NaN(Nu, Nsim);
uprev = us;

randn('state', 0);
v = 0*randn(Ny, Nsim + 1);
w = 0*randn(Nw, Nsim);

xhat = NaN(Nx, Nsim);
dhat = NaN(Nd, Nsim);

xtarg = NaN(Nx, Nsim);
utarg = NaN(Nu, Nsim);

xprior = repmat(xs, 1, Nmhe);
dprior = zeros(Nd, Nmhe);

% Disturbance and setpoint.
Nstart = 9;
p = [repmat(ps, 1, Nstart), repmat(1.1*ps, 1, Nsim - Nstart)];
ysp = NaN(Ny, Nsim + 1);
ysp(CVs,:) = repmat(ys(CVs), 1, Nsim + 1);

% Start loop.
for i = 1:(Nsim + 1)
    fprintf('(%3d) ', i);
    
    % Take measurement and run mhe.
    y(:,i) = C*x(:,i) + v(:,i);
    mhe.newmeasurement(y(:,i), uprev);
    mhe.par.xbar = xprior(:,1);
    mhe.par.dbar = dprior(:,1);
    mhe.solve();
    fprintf('Estimator: %s, ', mhe.status);
    if ~isequal(mhe.status, 'Solve_Succeeded')
        fprintf('\n');
        warning('mhe failed at time %d!', i);
        break
    end
    xhat(:,i) = mhe.var.x(:,end);
    dhat(:,i) = mhe.var.d(:,end);
    mhe.saveguess();
    xprior = [xprior(:,2:end), xhat(:,i)];
    dprior = [dprior(:,2:end), dhat(:,i)];
    
    % Stop if at last time.
    if i == Nsim + 1
        fprintf('Done\n');
        break
    end
    
    % Use steady-state target selector.
    sstarg.fixvar('y', 1, ysp(CVs,i), CVs);
    sstarg.par.d = dhat(:,i);
    sstarg.solve();
    fprintf('Target: %s, ', sstarg.status);
    if ~isequal(sstarg.status, 'Solve_Succeeded')
        fprintf('\n');
        warning('sstarg failed at time %d!', i);
        break
    end
    xtarg(:,i) = sstarg.var.x; 
    utarg(:,i) = sstarg.var.u;
    
    % Apply control law.
    controller.fixvar('x', 1, xhat(:,i));
    controller.par.xsp = xtarg(:,i);
    controller.par.usp = utarg(:,i);
    controller.par.d = dhat(:,i);
    if nonlinearMPC
        % Simulate to get a feasible initial guess.
        xguess = controller.guess.x;
        uguess = controller.guess.u;
        for t = 1:Ncontroller
            xguess(:,t + 1) = full(f(xguess(:,t), uguess(:,t), ...
                                     Bd, Bpd, dhat(:,i)));
        end
        controller.guess.x = xguess;
    end
    controller.solve();
    fprintf('Controller: %s, ', controller.status);
    if ~isequal(controller.status, 'Solve_Succeeded')
        fprintf('\n');
        warning('controller failed at time %d', i);
        break
    end
    
    u(:,i) = controller.var.u(:,1);
    controller.saveguess();
    uprev = u(:,i); % Save previous u.
    
    % Evolve plant.
    x(:,i + 1) = full(cstrsim(x(:,i), u(:,i), p(:,i))) + w(:,i);
    
    fprintf('\n');
end

% Make a plot.
if nonlinearMPC
    mpctype = 'Nonlinear';
else
    mpctype = 'Linear';
end
plottitle = sprintf('%s MPC: %s Disturbance Model', mpctype, disturbancemodel);
style = struct('fig', figure());
mpc.mpcplot('x', x, 'u', u, 'xnames', {'c', 'T', 'h'}, ...
            'unames', {'T_c', 'F'}, 'title', plottitle, 'legend', 'Actual', ...
            '**', style);
mpc.mpcplot('x', xhat, 'u', NaN(Nu, Nsim), 'color', 'b', 'linestyle', '--', ...
            'legend', 'Estimated', '**', style);
mpc.mpcplot('x', ysp, 'u', NaN(Nu, Nsim), 'color', 'r', 'linestyle', ':', ...
            'legend', 'Setpoint', '**', style);

