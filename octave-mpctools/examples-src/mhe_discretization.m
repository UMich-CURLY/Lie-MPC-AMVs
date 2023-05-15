% Compares MHE results using different discretizations for a moderately stiff
% nonisothermal CSTR system (From Rawlings and Ekerdt).
%
% Because of the difficult nonlinearities (in particular, reactor ignition), you
% need a pretty good initial guess. Keeping `use_perfect_guess` as true provides
% the actual states as a guess to MHE. Note, however, that because the various
% model discretizations are not exact, the guesses are not exactly "perfect,"
% and the solver still has to find a nearby solution that is feasible to the
% particular discretized model.
mpc = import_mpctools();
use_perfect_guess = true();

% Problem size.
Nx = 3;
Ny = 2;
Nw = Nx;
Nv = Ny;
Delta = 1;

% System model.
pars = struct();
pars.Tf = 298; % K
pars.Ta = 298; % K
pars.Tm = 298; % K
pars.Cphat = 4; % kJ/kg/K
pars.cAf = 2; % kmol/m^3
pars.cBf = 0;
pars.km = 0.004; % min^-1
pars.E = 1.5E4; % K
pars.rhof = 1e3; % kg/m^3
pars.dHr = -1.75e5; % kJ/kmol
pars.U0 = 340;
pars.tau = 72.3; % min
pars.nu = [-1; 1];
pars.C = [1 1 0; 0 0 1];

function dxdt = cstr_ode_pars(x, pars)
    cA = x(1); 
    cB = x(2);
    T = x(3);
    
    k = pars.km*exp(-pars.E*(1/T - 1/pars.Tm));
    Q_rxn = -pars.dHr*k*cA/(pars.rhof*pars.Cphat);
    Q_flow = (pars.Tf - T)/pars.tau;
    Q_cooling = pars.U0*(pars.Ta - T)/(pars.rhof*pars.Cphat);

    dxdt = [(pars.cAf - cA)/pars.tau + pars.nu(1)*k*cA;
            (pars.cBf - cB)/pars.tau + pars.nu(2)*k*cA;
            Q_rxn + Q_flow + Q_cooling];
end%function
cstr_ode = @(x) cstr_ode_pars(x, pars);
simulator = mpc.getCasadiIntegrator(cstr_ode, Delta, ...
                                    [Nx], {'x'}, {'simulator'});

function y = cstr_measurement_pars(x, pars)
    y = pars.C*x;
end%function
cstr_measurement = @(x) cstr_measurement_pars(x, pars);
h = mpc.getCasadiFunc(cstr_measurement, [Nx], {'x'}, {'h'});

% Variable bounds.
lb = struct('x', [0; 0; 280]);
ub = struct('x', [2.1; 2.1; 390]);

% MHE weighting.
Nt = 250;
P = diag(([2.5, 2.5, 250]/Nt).^2);
Q = diag(([0.005, 0.005, 1]*Delta).^2);
R = diag(([0.25, 0.2]*Delta).^2);

[Pinv, Phalf] = mpctools.spdinv(P);
[Qinv, Qhalf] = mpctools.spdinv(Q);
[Rinv, Rhalf] = mpctools.spdinv(R);

l = mpc.getCasadiFunc(@(w, v) w'*Qinv*w + v'*Rinv*v, ...
                      [Nw, Nv], {'w', 'v'}, {'l'});
lx = mpc.getCasadiFunc(@(x, x0bar) (x - x0bar)'*Pinv*(x - x0bar), ...
                       [Nx, Nx], {'x', 'x0bar'}, {'lx'});

% Simulate to get data.
randn('state', 1);
wsim = Qhalf*randn(Nw, Nt);
vsim = Rhalf*randn(Nv, Nt + 1);

xsim = NaN(Nx, Nt + 1);
xsim(:,1) = [pars.cAf; 0.001; pars.Tf];

ysim = NaN(Ny, Nt + 1);
yclean = NaN(Ny, Nt + 1);
for t = 1:(Nt + 1)
    yclean(:,t) = full(h(xsim(:,t)));
    ysim(:,t) = yclean(:,t) + vsim(:,t);
    if t <= Nt
        xsim(:,t + 1) = full(simulator(xsim(:,t))) + wsim(:,t);
    end
end

% Define explicit Euler discretization.
function x = euler(f, Delta, M, x)
    % Simulates explicit Euler method.
    narginchk(4, 4);
    Delta = Delta/M;
    for i = 1:M
        x = x + Delta*f(x);
    end
end%function

% Test different f discretizations.
N = struct('x', Nx, 'y', Ny, 'w', Nw, 'v', Nv, 't', Nt);
discretizations = {
    'euler', 1;
    'euler', 3;
    'rk4', 1;
    'rk4', 2;
    'colloc', 2;
    'colloc', 4;
    'integrator', 1;
};
data = struct();
for i = 1:size(discretizations, 1)
    name = sprintf('%s_%d', discretizations{i,:});
    fprintf('Using %s: ', name);
    par = struct('y', ysim, 'x0bar', xsim(:,1));
    if use_perfect_guess
        guess = struct('x', xsim, 'y', yclean, ...
                       'w', wsim, 'v', vsim);
    else
        guess = struct('x', 0.5*(lb.x + ub.x));
        guess.y = cstr_measurement(guess.x);
    end
    kwargs = struct('h', h, 'l', l, 'lx', lx, 'N', N, 'wadditive', true());
    
    % Define discretized f function and build estimator.
    switch discretizations{i,1}
    case 'euler'
        f = @(x) euler(cstr_ode, Delta, discretizations{i,2}, x);
        f = mpc.getCasadiFunc(f, [Nx], {'x'}, {['f_', name]});
    case 'rk4'
        f = mpc.getCasadiFunc(cstr_ode, [Nx], {'x'}, {['f_', name]}, ...
                              'rk4', true(), 'Delta', Delta, ...
                              'M', discretizations{i,2});
    case 'colloc'
        f = mpc.getCasadiFunc(cstr_ode, [Nx], {'x'}, {['f_', name]});
        kwargs.Delta = Delta;
        kwargs.N.c = discretizations{i,2};
    case 'integrator'
        f = simulator; % Use actual simulator.
        kwargs.casaditype = 'MX';
    otherwise
        error('Invalid discretization choice!');
    end
    buildtime = tic();
    estimator = mpc.nmhe('f', f, 'par', par, 'lb', lb, 'ub', ub, ...
                         'guess', guess, '**', kwargs);
    buildtime = toc(buildtime);
    fprintf('(%.2f s build) ', buildtime);
    
    % Simulate system.
    solvetime = tic();
    estimator.solve();
    solvetime = toc(solvetime);
    fprintf('(%.2f s solve) %s\n', solvetime, estimator.status);
    data.(name) = estimator.var;
    data.(name).y = par.y - estimator.var.v;
    
    % Calculate RMSE of concentration measurements.
    Cerr = data.(name).x(1:2,:) - xsim(1:2,:);
    rmse = sqrt(mean(sum(Cerr.^2, 1)));
    fprintf('    RMSE of cA, cB estimates: %g\n', rmse);
end


% Make a plot.
fig = figure();
t = Delta*(0:Nt);
kwargs = struct('fig', figure(), 'xnames', {{'c_A', 'c_B', 'T'}}, ...
            'unames', {{'c_A + c_B', 'T'}});
mpc.mpcplot(xsim, ysim, t, 'legend', 'actual', 'color', 'k', ...
            'marker', 'x', '**', kwargs);
keys = fieldnames(data);
colors = jet(length(keys));
markers = 'osvd^v<>ph';
for i = 1:length(keys)
    key = keys{i};
    mpc.mpcplot(data.(key).x, data.(key).y, t, ...
                'legend', key, 'color', colors(i,:), ...
                'marker', markers(i), '**', kwargs);
end

