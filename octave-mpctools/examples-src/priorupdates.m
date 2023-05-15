% Comparison of different prior update strategies for a nonlinear system.
mpc = import_mpctools();
use_collocation = true();

Nx = 2;
Ny = 1;
Nt = 10; % MHE horizon.

k1 = 0.16;
k2 = 0.0064;

Delta = 0.1;

% Define system model.
ode = @(x) [-2*k1*x(1)^2 + 2*k2*x(2); ...
            k1*x(1)^2 - k2*x(2)];
measurement = @(x) x(1) + x(2);

sim = mpc.getCasadiIntegrator(ode, Delta, [Nx], {'x'}, 'funcname', 'sim');
if use_collocation
    f = mpc.getCasadiFunc(ode, [Nx], {'x'}, 'funcname', 'f');
else
    f = mpc.getCasadiFunc(ode, [Nx], {'x'}, 'funcname', 'f', ...
                          'rk4', true(), 'Delta', Delta);
end
h = mpc.getCasadiFunc(measurement, [Nx], {'x'}, 'funcname', 'h');

% Simulate the system.
Nsim = 100;
xsim = NaN(Nx, Nsim + 1);
xsim(:,1) = [3; 1];
ysim = NaN(Ny, Nsim + 1);
for i = 1:(Nsim + 1)
    ysim(:,i) = measurement(xsim(:,i));
    if i <= Nsim
        xsim(:,i + 1) = full(sim(xsim(:,i)));
    end
end

% Choose MHE objective function.
Q = diag([0.001, 0.01].^2);
R = diag([0.1].^2);
P0 = diag([1, 1].^2);

Qinv = mpc.spdinv(Q);
Rinv = mpc.spdinv(R);
P0inv = mpc.spdinv(P0);

l = mpc.getCasadiFunc(@(w, v) w'*Qinv*w + v'*Rinv*v, [Nx, Ny], {'w', 'v'}, ...
                      'funcname', 'l');

par = struct('x0bar', [0.1; 4.5], 'Pinv', P0inv); % Initial prior.

% Choose other parameters.
lb = struct('x', zeros(Nx, 1));
N = struct('x', Nx, 'y', Ny, 't', Nt);
if use_collocation
    N.c = 2;
end
kwargs = struct('f', f, 'h', h, 'l', l, 'N', N, 'lb', lb, 'par', par, ...
                'y', ysim(:,1:(N.t + 1)), 'wadditive', true());
if use_collocation
    kwargs.Delta = Delta;
end

% Simulate different prior updates.
updates = {'filtering', 'smoothing'};
mhes = struct();
xmhes = struct('actual', xsim);
ymhes = struct('actual', ysim);
for k = 1:length(updates)
    key = updates{k};
    fprintf('*** Simulating MHE with %s update.\n', key);
    mhe = mpc.nmhe('priorupdate', key, '**', kwargs);
    mhe.fix_truncated_x = true();
    
    xmhe = NaN(Nx, Nsim + 1);
    ymhe = NaN(Ny, Nsim + 1);
    for i = 0:Nsim
        if i <= Nt
            mhe.truncatehorizon(i);
        else
            mhe.newmeasurement(ysim(:,i + 1));
        end
        
        mhe.solve();
        fprintf('i = %d: %s\n', i, mhe.status);
        if ~isequal(mhe.status, 'Solve_Succeeded')
            warning('Solver failed!');
            break
        end
        mhe.saveestimate();
        xmhe(:,i + 1) = mhe.history(1).xhat(:,end);
        ymhe(:,i + 1) = mhe.history(1).yhat(:,end);
    end
    
    mhes.(key) = mhe;
    xmhes.(key) = xmhe;
    ymhes.(key) = ymhe;
end

% Make a plot.
keys = {'actual', 'filtering', 'smoothing'};
colors = {'k', 'r', 'g'};
markers = {'x', 's', 'o'};
plotstyle = struct('t', (0:Nsim)*Delta, 'fig', figure(), ...
                   'xnames', {{'p_a', 'p_b', 'p_{tot}'}}, ...
                   'unames', {{'log_{10}|e_a|', 'log_{10}|e_b|'}});
for i = 1:length(keys)
    k = keys{i};
    x = [xmhes.(k); ymhes.(k)];
    e = log10(abs(xmhes.(k) - xmhes.actual));
    mpc.mpcplot(x, e, 'legend', k, 'color', colors{i}, ...
                'marker', markers{i}, '**', plotstyle);
end

