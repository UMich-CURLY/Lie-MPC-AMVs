% Test QP export on a linear problem.
mpc = import_mpctools();

% Define parameters and system model.
T = 1;
Delta = 0.01;
Nx = 1;
Nu = 1;
Nt = 100;

Ac = -1;
Bc = 10;
model = mpc.c2d(Delta, Ac, Bc);
f = @(x, u) model.A*x + model.B*u; % Model.
fcasadi = mpc.getCasadiFunc(f, [Nx, Nu], {'x', 'u'}, {'F'});

% Pick time-varying setpoint.
T = 0:0.25:ceil(Nt*Delta/T);
Xsp = zeros(size(T));
Xsp(2:4:end) = 1;
Xsp(4:4:end) = -1;

t = (0:Nt)*Delta;
xsp = interp1(T, Xsp, t);
usp = zeros(Nu, Nt); % This is an inconsistent setpoint.

% Pick time-varying weights. Time runs along third dimension.
lam = reshape(1 + cos(2*pi()*t(1:(end - 1))/t(end)), 1, 1, Nt)/2;
Q = bsxfun(@times, eye(Nx), lam);
R = bsxfun(@times, eye(Nu), 1 - lam);

% Define stage cost.
l = @(x, u, xsp, usp, Q, R) (x - xsp)'*Q*(x - xsp) + (u - usp)'*R*(u - usp);
lcasadi = mpc.getCasadiFunc(l, {Nx, Nu, Nx, Nu, [Nx, Nx], [Nu, Nu]}, ...
                            {'x', 'u', 'xsp', 'usp', 'Q', 'R'}, {'l'});

% Build bounds, parameters, and N.
lb = struct('u', -ones(Nu, Nt));
ub = struct('u', ones(Nu, Nt));
par = struct('xsp', xsp, 'usp', usp, 'Q', Q, 'R', R);

N = struct('x', Nx, 'u', Nu, 't', Nt);

x0 = -2;

% Build solver and optimize.
solver = mpc.nmpc('f', fcasadi, 'l', lcasadi, 'N', N, 'lb', lb, 'ub', ub, ...
                  'x0', x0, 'par', par, 'verbosity', 0);

% Solve using exported solvers and native solve.
sols = struct();
if exist('gurobi')
    prob = solver.getQP('gurobi');
    sol = gurobi(prob);
    sols.gurobi = solver.xvec2struct(sol.x);
else
    fprintf('Skipping gurobi.\n');
end
if exist('gurobi_octave')
    prob = solver.getQP('gurobi');
    sol = gurobi_octave(prob);
    sols.gurobioctave = solver.xvec2struct(sol.x);
else
    fprintf('Skipping gurobi_octave.\n');
end
if exist('qp')
    args = solver.getQP('qp');
    [z, ~, status] = qp(args{:});
    if status.info ~= 0
        warning('qp returned status %d!', status.info);
    end
    sols.qp = solver.xvec2struct(z);
else
    fprintf('Skipping qp.\n');
end

solver.solve();
sols.ipopt = solver.var;

% Make a plot.
options = struct('fig', figure(), 'umarker', 'yes');
keys = fieldnames(sols);
colors = jet(length(keys));
markers = 'os^dv568<>';
for i = 1:length(keys)
    k = keys{i};
    x = sols.(k).x;
    u = sols.(k).u;
    mpc.mpcplot([x; x - sols.ipopt.x], [u; u - sols.ipopt.u], ...
                'color', colors(i,:), 'legend', k, ...
                'xnames', {'x', '\Delta x'}, 'unames', {'u', '\Delta u'}, ...
                'marker', markers(i), '**', options);
end

% Check solutions.
checks = setdiff(fieldnames(sols), 'ipopt');
for i = 1:length(checks)
    c = checks{i};
    assert(sols.ipopt.x, sols.(c).x, 1e-6);
    assert(sols.ipopt.u, sols.(c).u, 1e-6);
end

