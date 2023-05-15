% Test MPCTools on a simple linear MPC problem.
mpc = import_mpctools();

% Define system.
A = 1;
B = 1;
Q = 1;
R = 1;
x0 = 1;

pkg('load','control');
[K, P] = dlqr(A, B, Q, R);

f = @(x, u) A*x + B*u;
l = @(x, u) x'*Q*x + u'*R*u;
Vf = @(x) x'*P*x;

N = struct('x', 1, 'u', 1, 't', 10);

% Create Casadi functions.
fcasadi = mpc.getCasadiFunc(f, [N.x, N.u], {'x', 'u'}, {'f'});
lcasadi = mpc.getCasadiFunc(l, [N.x, N.u], {'x', 'u'}, {'l'});
Vfcasadi = mpc.getCasadiFunc(Vf, [N.x], {'x'}, {'Vf'});

% Solve with IPOPT.
solver = mpc.nmpc('f', fcasadi, 'N', N, 'x0', x0, 'l', lcasadi, ...
                  'Vf', Vfcasadi, 'verbosity', 5);
solver.solve();
sol = solver.var;

% Double-check solution.
check = struct('x', NaN(N.x, N.t + 1), 'u', NaN(N.u, N.t));
check.x(:,1) = x0;
for i = 1:N.t
    check.u(:,i) = -K*check.x(:,i);
    check.x(:,i + 1) = A*check.x(:,i) + B*check.u(:,i);
end

fprintf('MPC cost: %g\n', solver.obj);
fprintf('LQR cost: %g\n', x0'*P*x0);

% Compare.
kwargs = struct('fig', figure(), 'legendloc', 'SouthEast');
mpc.mpcplot(check.x, check.u, 'legend', 'LQR', 'color', 'r', 'marker', 'x', ...
            '**', kwargs);
mpc.mpcplot(sol.x, sol.u, 'legend', 'MPC', 'color', 'g', 'marker', 'o', ...
            '**', kwargs);
