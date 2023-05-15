% Applies offset-free linear MPC to the nonlinear CSTR.
% See Pannocchia and Rawlings, AIChE J, 2002.
pkg('load', 'control');
mpc = import_mpctools();

% Parameters and sizes for the nonlinear system
Delta = 1;
Nx = 3;
Nu = 2;
Ny = Nx;
Np = 1;
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

ode_casadi = mpc.getCasadiFunc(ode, [Nx, Nu, Np], ...
                               {'x', 'u', 'p'}, {'ode'});
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

% Simulate a few steps so we actually get dx/dt = 0 at steady state.
for i = 1:10
    xs = full(cstrsim(xs, us, ps));
end
cs = xs(1);
Ts = xs(2);
hs = xs(3);

% Get linearized model and linear controller.
model = mpc.getLinearizedModel(ode_casadi, {xs, us, ps}, ...
                               {'A', 'B', 'Bp'}, Delta);
A = model.A;
B = model.B;
C = eye(Ny);
Bp = model.Bp;
Q = diag(1./xs.^2);
R = diag(1./us.^2);
[K, Pi] = dlqr(A, B, Q, R);
K = -K;

% Pick whether to use good disturbance model.
useGoodDisturbanceModel = true();

if useGoodDisturbanceModel
    % disturbance model 6; no offset
    Nd = 3;
    Bd = zeros(Nx, Nd);
    Bd(:,3) = B(:,2);
    Cd = [1 0 0; 0 0 0; 0 1 0];
else
    % disturbance model with offset
    Nd = 2;
    Bd = zeros(Nx, Nd);
    Cd = [1 0; 0 0; 0 1];
end

% Set up state estimator; use KF
Qw = zeros(Nx + Nd);
Qw(1:Nx,1:Nx) = small*eye(Nx);
Qw(Nx+1:end,Nx+1:end) = small*eye(Nd);
Qw(end,end) = 1.0;
Rv = small*diag(xs.^2);

kf = mpc.KalmanFilter('A', A, 'B', B, 'C', C, 'Bd', Bd, 'Cd', Cd, ...
                      'Qw', Qw, 'Rv', Rv, 'contvars', [1, 3]);

% Closed-loop simulation.
Nsim = 50;
x = NaN(Nx, Nsim + 1);
x(:, 1) = 0;
y = NaN(Ny, Nsim + 1);
u = NaN(Nu, Nsim);

v = zeros(Ny, Nsim + 1);

xhatm = NaN(Nx, Nsim + 1);
xhatm(:,1) = 0;
dhatm = NaN(Nd, Nsim + 1);
dhatm(:,1) = 0;

xhat = NaN(Nx, Nsim);
dhat = NaN(Nd, Nsim);

xtarg = zeros(Nx, Nsim);
utarg = zeros(Nu, Nsim);

% Disturbance and setpoint.
p = zeros(Np, Nsim);
p(10:end) = 0.1*F0s;
ysp = zeros(Ny, Nsim);

% Start loop.
for i = 1:(Nsim + 1)
    fprintf('Step %i of %d\n', i - 1, Nsim);
    
    % Take measurement.
    y(:,i) = C*x(:,i) + v(:,i);
    
    % Advance state measurement.
    [xhat(:,i), dhat(:,i)] = kf.filter(y(:,i), xhatm(:,i), dhatm(:,i));
    
    % Stop if at last time.
    if i == Nsim + 1
        break
    end
    
    % Use steady-state target selector.
    [xtarg(:,i), utarg(:,i)] = kf.target(ysp(:,i));
    
    % Apply control law.
    u(:,i) = K*(xhat(:,i) - xtarg(:,i)) + utarg(:,i);
    
    % Evolve plant. Our variables are deviation but cstrsim needs positional.
    x(:,i + 1) = full(cstrsim(x(:,i) + xs, u(:,i) + us, p(:,i) + ps)) - xs;
    
    % Advance state estimates
    [xhatm(:,i + 1), dhatm(:,i + 1)] = kf.predict(u(:,i), xhat(:,i), dhat(:,i));
end
u(:,end) = u(:,end-1); % Repeat for stair plot.

% Make a plot.
mpc.mpcplot('x', x, 'u', u, 'xnames', {'c', 'T', 'h'}, 'unames', {'T_c', 'F'});

% <--
% Save data.
save('-v7', 'cstr.mat', 'x', 'u');
% -->

