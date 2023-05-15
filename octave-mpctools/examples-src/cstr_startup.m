% Startup of a nonlinear CSTR using LMPC and NMPC.
pkg('load', 'control');
mpc = import_mpctools();

% Parameters and sizes for the nonlinear system
Delta = 0.25;
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
ode_rk4_casadi = mpc.getCasadiFunc(ode, [Nx, Nu, Np], {'x', 'u', 'p'}, ...
                                   {'ode_rk4'}, 'rk4', true(), 'Delta', Delta);
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
[K, Pi]  = dlqr(A, B, Q, R);
K = -K;

% Define Casadi functions.
Flinear = @(x,u,p) model.A*(x - xs) + model.B*(u - us) + model.Bp*(p - ps) + xs;
Flinear = mpc.getCasadiFunc(Flinear, [Nx, Nu, Np], {'x', 'u', 'p'}, {'F'});
Fnonlinear = ode_rk4_casadi;

l = @(x, u, xsp, usp) (x - xsp)'*Q*(x - xsp) + (u - usp)'*R*(u - usp);
l = mpc.getCasadiFunc(l, [Nx, Nu, Nx, Nu], {'x', 'u', 'xsp', 'usp'}, {'l'});

h = mpc.getCasadiFunc(@(x,p) x, [Nx, Np], {'x', 'p'}, {'h'});

Vf = mpc.getCasadiFunc(@(x, xsp) (x - xsp)'*Pi*(x - xsp), [Nx, Nx], ...
                       {'x', 'xsp'}, {'Vf'});

% Choose simulation time and preallocate everybody.
Nsim = 100;
x0 = [0.05*cs; 0.75*Ts; 0.5*hs];
xcl = struct();
ucl = struct();

% See what happens if we just try to use the steady-state input.
xcl.uncont = NaN(Nx, Nsim + 1);
xcl.uncont(:,1) = x0;
ucl.uncont = repmat(us, 1, Nsim);
for t = 1:Nsim
    xcl.uncont(:,t + 1) = full(cstrsim(xcl.uncont(:,t), ucl.uncont(:,t), ps));
end

% Build solvers for linear and nonlinear models.
Nt = 15;

umax = [0.05*Tcs; 0.15*Fs];
Dumax = 0.2*umax;
ulb = us - umax;
uub = us + umax;

lb = struct();
lb.u = repmat(ulb, 1, Nt);
lb.Du = repmat(-Dumax, 1, Nt);

ub = struct();
ub.u = repmat(uub, 1, Nt);
ub.Du = repmat(Dumax, 1, Nt);

par = struct();
par.xsp = repmat(xs, 1, Nt + 1);
par.usp = repmat(us, 1, Nt);
par.p = repmat(ps, 1, Nt + 1);
par.uprev = us;

N = struct('x', Nx, 'u', Nu, 'p', Np, 't', Nt, 'y', Ny);

kwargs = struct();
kwargs.N = N;
kwargs.l = l;
kwargs.Vf = Vf;
kwargs.lb = lb;
kwargs.ub = ub;
kwargs.par = par;

solvers = struct();
solvers.lmpc = mpc.nmpc('f', Flinear, '**', kwargs);
solvers.nmpc = mpc.nmpc('f', Fnonlinear, '**', kwargs);

% Build steady-state target finder.
contvars = [1; 3]; % Temperature and height are controlled.
kwargs = struct();
kwargs.N = N;
kwargs.lb = struct('u', us - umax, 'x', zeros(Nx, 1));
kwargs.ub = struct('u', us + umax);
kwargs.h = h;
kwargs.par = struct('p', ps);

sstargs = struct();
sstargs.lmpc = mpc.sstarg('f', Flinear, '**', kwargs);
sstargs.nmpc = mpc.sstarg('f', Fnonlinear, '**', kwargs);

% Simulate closed-loop control.
solvernames = fieldnames(solvers);
for i = 1:length(solvernames)
    s = solvernames{i};
    solver = solvers.(s);
    sstarg = sstargs.(s);
    
    solver.saveguess(struct('x', par.xsp, 'u', par.usp));
    solver.par.uprev = us;
    
    xcl.(s) = NaN(Nx, Nsim + 1);
    xcl.(s)(:,1) = x0;
    
    ucl.(s) = NaN(Nu, Nsim);
    
    ysp = repmat(xs, 1, Nsim + 1);
    changesp = round(Nsim/3):round(2*Nsim/3);
    newsp = repmat(xs.*[0.85; 0.75; 1.15], 1, Nsim + 1);
    ysp(contvars,changesp) = newsp(contvars,changesp);
    
    for t = 1:Nsim
        if t == 1 || ~isequal(ysp(:,t), ysp(:,t - 1))
            % Use steady-state target finder.
            sstarg.fixvar('y', 1, ysp(contvars,t), contvars);
            sstart.guess.y = ysp(:,t);
            sstarg.guess.x = ysp(:,t);
            sstarg.guess.u = us;
            sstarg.solve();
            fprintf('%10s %3d: %s\n', 'sstarg', t, sstarg.status);
            if ~isequal(sstarg.status, 'Solve_Succeeded')
                fprintf('*** sstarg failed!\n');
                break
            end
            
            % Set setpoints.
            solver.par.xsp = repmat(sstarg.var.x, 1, Nt + 1);
            solver.par.usp = repmat(sstarg.var.u, 1, Nt);
        end
        
        % Set initial condition and solve.
        solver.fixvar('x', 1, xcl.(s)(:,t));
        solver.solve();
        fprintf('%10s %3d: %s\n', s, t, solver.status);
        if ~isequal(solver.status, 'Solve_Succeeded')
            fprintf('*** Solver failed!\n');
            break
        end
        solver.saveguess()
        
        % Simulate system.
        ucl.(s)(:,t) = solver.var.u(:,1);
        xcl.(s)(:,t + 1) = full(cstrsim(xcl.(s)(:,t), ucl.(s)(:,t), ps));
        
        % Update previous u.
        solver.par.uprev = ucl.(s)(:,t);
    end
end

% Make a plot.
fig = figure();
xsp = NaN(Nx, Nsim);
t = Delta*(0:Nsim);
xsp(contvars,:) = ysp(contvars,1:end-1);
plotargs = struct('fig', fig, 't', t);
mpc.mpcplot(xcl.uncont, ucl.uncont, 'xsp', xsp, 'color', 'r', ...
            'spcolor', 'k', 'legend', 'Uncontrolled', '**', plotargs);
mpc.mpcplot(xcl.lmpc, ucl.lmpc, 'legend', 'LMPC', 'color', 'b', '**', plotargs);
mpc.mpcplot(xcl.nmpc, ucl.nmpc, 'legend', 'NMPC', 'color', 'g', ...
            'xnames', {'c', 'T', 'h'}, 'unames', {'T_c', 'F'}, '**', plotargs);

% <--
% Save data for plot.
save('-v7', 'cstr_startup.mat', 'xcl', 'ucl', 'xsp', 't', 'contvars', ...
     'ulb', 'uub');
% -->

