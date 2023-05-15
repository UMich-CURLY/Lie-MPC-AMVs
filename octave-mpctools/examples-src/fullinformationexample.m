% Full-information estimation on a nonlinear batch reactor.
mpc = import_mpctools();
rand('seed', 0);

% Parmeters.
Delta = 0.25;
Nsim = 20;

Nx = 3;
Nu = 1; % There isn't actually an input.
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

Pinv = mpctools.spdinv(P);
Qinv = mpctools.spdinv(Q);
Rinv = mpctools.spdinv(R);

xhat0 = [1; 0; 4];
x0 = [0.5; 0.05; 0.0];

% Model parameters.
pars = struct();
pars.k1 = 0.5;
pars.km1 = 0.05;
pars.k2 = 0.2;
pars.km2 = 0.01;
pars.RT = 32.84;

function dxdt = odefunc(x, u, w, pars)
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
ode = @(x, u, w) odefunc(x, u, w, pars);
model = mpc.getCasadiIntegrator(ode, Delta, [Nx, Nu, Nw], ...
                                {'x', 'u', 'w'}, {'cstr'});
F = mpc.getCasadiFunc(ode, [Nx, Nu, Nw], {'x', 'u', 'w'}, {'f'}, ...
                      'rk4', true(), 'Delta', Delta, 'M', 4);

function y = measfunc(x, pars)
    % Measurement function.
    y = pars.RT*(x(1) + x(2) + x(3));
end%function
meas = @(x) measfunc(x, pars);
H = mpc.getCasadiFunc(meas, [Nx], {'x'}, {'H'});

% Pick stage costs.
lfunc = @(w, v) w'*Qinv*w + v'*Rinv*v;
lxfunc = @(x, x0bar, Pinv) (x - x0bar)'*Pinv*(x - x0bar);

l = mpc.getCasadiFunc(lfunc, [Nw, Nv], {'w', 'v'}, {'l'});
lx = mpc.getCasadiFunc(lxfunc, {Nx, Nx, [Nx, Nx]}, ...
                       {'x', 'x0bar', 'Pinv'}, {'lx'});

% First, simulate the system to get data.
w = sig_w*randn(Nw, Nsim);
v = sig_v*randn(Nv, Nsim + 1);

xsim = NaN(Nx, Nsim + 1);
xsim(:,1) = x0;

usim = zeros(Nu, Nsim);

ysim = NaN(Ny, Nsim + 1);
yclean = NaN(Ny, Nsim + 1);

for t = 1:(Nsim + 1)
    yclean(:,t) = full(H(xsim(:,t)));
    ysim(:,t) = yclean(:,t) + v(:,t);
    
    if t <= Nsim
        xsim(:,t + 1) = full(model(xsim(:,t), usim(:,t), w(:,t)));
    end
end

% Now build MHE solver.
par = struct('Pinv', Pinv);

N = struct('x', Nx, 'u', Nu, 'w', Nw, 'y', Ny);

x0bar = xhat0;

lb = struct(); % Will fill later.

nmheargs = struct('f', F, 'l', l, 'h', H, 'N', N, 'lx', lx, ...
                  'x0bar', x0bar, 'par', par, ...
                  'lb', struct(), 'guess', struct());

% Initial step is odd because there is no model. Call nlfilter.
xhat = NaN(Nx, Nsim + 1, Nsim + 1);
xplot = NaN(Nx, Nsim + 1);
yhat = NaN(Ny, Nsim + 1, Nsim + 1);
yplot = NaN(Ny, Nsim + 1);
[xhat(:,1,1), status] = mpc.nlfilter(H, x0bar, ysim(:,1), Rinv, Pinv, ...
                                     'xlb', zeros(Nx, 1), 'xub', inf(Nx, 1));
yhat(:,1,1) = full(H(xhat(:,1,1)));
fprintf('Step 0: %s\n', status);
xplot(:,1) = xhat(:,1,1);
yplot(:,1) = yhat(:,1,1);

% For the remaining steps, we have to build a new solver object since the
% horizon is changing.
totaltime = tic();
for t = 1:Nsim
    % Adjust time-varying arguments.
    nmheargs.N.t = t;
    nmheargs.y = ysim(:,1:t + 1);
    nmheargs.u = usim(:,1:t);
    nmheargs.lb.x = zeros(Nx, t + 1);
    
    % Solve the optimization problem.
    buildtime = tic();
    solver = mpc.nmhe('**', nmheargs);
    buildtime = toc(buildtime);
    
    solvetime = tic();
    solver.solve();
    solvetime = toc(solvetime);
    
    fprintf('%3d (%5.3g s build, %5.3g s solve): %s\n', t, ...
            buildtime, solvetime, solver.status);
    
    % Save the optimal vales and adjust the guess.
    xhat(:, 1:t + 1, t + 1) = solver.var.x;
    yhat(:, 1:t + 1, t + 1) = solver.par.y - solver.var.v;
    xplot(:,t + 1) = xhat(:,t + 1,t + 1);
    yplot(:,t + 1) = yhat(:,t + 1,t + 1);
    
    guess = solver.var;
    guessfields = fieldnames(guess);
    for i = 1:length(guessfields)
        f = guessfields{i};
        guess.(f) = [guess.(f), guess.(f)(:,end)];
    end
    nmheargs.guess = guess;
end
totaltime = toc(totaltime);
fprintf('Simulation took %.5g s\n', totaltime);

% Make a plot.
colors = {'r', 'b', 'g'};
figure();

subplot(2, 1, 1);
hold('on');
for i = 1:Nx
    % Plot initial sequence, best current estimates, and actual data.
    plot((0:Nsim)*Delta, xplot(i,:), ['-', colors{i}]);
    plot((0:Nsim)*Delta, xsim(i,:), ['o', colors{i}]);
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


