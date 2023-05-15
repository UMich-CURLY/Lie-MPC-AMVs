% MHE on a nonlinear batch reactor. From Example 4.27 (Rawlings and Mayne, 2009).
mpc = import_mpctools();
randn('seed', 0);

% Parmeters.
Nt = 10;
Delta = 0.25;
Nsim = 80;

Nx = 3;
Ny = 1;
Nw = Nx;
Nv = Ny;

% Random variable standard deviations.
sig_v = 0.25; % Measurement noise.
sig_w = 0.001; % State noise.
sig_p = 0.5; % Prior.

P = sig_p.^2*eye(Nx);
Q = sig_w.^2*eye(Nw);
R = sig_v.^2*eye(Ny);

xhat0 = [1; 0; 4];
x0 = [0.5; 0.05; 0.0];

% Model parameters.
pars = struct();
pars.k1 = 0.5;
pars.km1 = 0.05;
pars.k2 = 0.2;
pars.km2 = 0.01;
pars.RT = 32.84;

function dxdt = odefunc(x, w, pars)
    % Continuous-time model.
    cA = x(1);
    cB = x(2);
    cC = x(3);
    rate1 = pars.k1*cA - pars.km1*cB*cC;
    rate2 = pars.k2*cB.^2 - pars.km2*cC;
    dxdt = [
        -rate1 + w(1);
        rate1 - 2*rate2 + w(2);
        rate1 + rate2 + w(3);
    ];
end%function
ode = @(x, w) odefunc(x, w, pars);
model = mpc.getCasadiIntegrator(ode, Delta, [Nx, Nw], ...
                                {'x', 'w'}, {'cstr'});
F = mpc.getCasadiFunc(ode, [Nx, Nw], {'x', 'w'}, {'f'}, ...
                      'rk4', true(), 'Delta', Delta, 'M', 4);

function y = measfunc(x, pars)
    % Measurement function.
    y = pars.RT*(x(1) + x(2) + x(3));
end%function
meas = @(x) measfunc(x, pars);
H = mpc.getCasadiFunc(meas, [Nx], {'x'}, {'H'});

% Pick stage costs.
lfunc = @(w, v) (w'*w)/sig_w.^2 + v'*v/sig_v.^2;

l = mpc.getCasadiFunc(lfunc, [Nw, Nv], {'w', 'v'}, {'l'});

% First, simulate the system to get data.
w = sig_w*randn(Nw, Nsim);
v = sig_v*randn(Nv, Nsim + 1);

xsim = NaN(Nx, Nsim + 1);
xsim(:,1) = x0;

ysim = NaN(Ny, Nsim + 1);
yclean = NaN(Ny, Nsim + 1);

for t = 1:(Nsim + 1)
    yclean(:,t) = full(H(xsim(:,t)));
    ysim(:,t) = yclean(:,t) + v(:,t);
    
    if t <= Nsim
        xsim(:,t + 1) = full(model(xsim(:,t), w(:,t)));
    end
end

% Now build MHE solver.
par = struct('Pinv', mpctools.spdinv(P), 'x0bar', xhat0);

N = struct('x', Nx, 'w', Nw, 'y', Ny, 't', Nt);

lb = struct('x', zeros(Nx, Nt + 1)); % Concentrations are nonnegative.
ub = struct('x', 10*ones(Nx, Nt + 1)); % Concentrations shouldn't be this high.

buildsolvertime = tic();
solver = mpc.nmhe('f', F, 'l', l, 'h', H, 'y', ysim(:,1:Nt + 1), 'N', N, ...
                  'lb', lb, 'ub', ub, 'par', par, 'Nhistory', Nsim, ...
                  'priorupdate', 'filtering');
buildsolvertime = toc(buildsolvertime);
fprintf('Building solver took %g s.\n', buildsolvertime);

% Loop.
xhat = NaN(Nx, Nt + 1, Nsim + 1);
xplot = NaN(Nx, Nsim + 1);
yhat = NaN(Ny, Nt + 1, Nsim + 1);
yplot = NaN(Ny, Nsim + 1);
tic();
for t = 1:(Nsim + 1)
    % Get new measurement or extend horizon.
    if t > Nt + 1
        solver.newmeasurement(ysim(:,t));
    else
        solver.truncatehorizon(t - 1);
    end
    
    % Solve MHE problem and save state estimate.
    solver.solve();
    fprintf('Step %d: %s\n', t, solver.status);
    if ~isequal(solver.status, 'Solve_Succeeded')
        warning('Solver failure at time %d!', t);
        break
    end
    solver.saveestimate(); % Stores current estimate to struct.
    xplot(:,t) = solver.history(1).xhat(:,end);
    yplot(:,t) = solver.history(1).yhat(:,end);
end
fprintf('Took %.4g s.\n', toc());

% Make a plot.
colors = {'r', 'b', 'g'};
figure();

subplot(2, 1, 1);
hold('on');
for i = 1:Nx
    % Plot initial sequence, best current estimates, and actual data.
    plot((0:Nsim)*Delta, xplot(i,:), ['-', colors{i}]);
    plot((0:Nsim)*Delta, squeeze(xsim(i,:)), ['o', colors{i}]);
end
ylabel('Concentration');
legend('Est. c_A', 'Act. c_A', 'Est. c_B', 'Act. c_B', 'Est. c_C', ...
       'Act. c_C', 'location', 'EastOutside');

subplot(2, 1, 2);
plot((0:Nsim)*Delta, yplot, '-k', (0:Nsim)*Delta, ysim, 'xk', ...
     (0:Nsim)*Delta, yclean, 'ok');
ylabel('Pressure');
xlabel('Time');
legend('Est. P', 'Act. P', 'Meas. P', 'location', 'EastOutside');

