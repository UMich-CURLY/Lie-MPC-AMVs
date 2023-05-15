% MHE applied to a linear system.
mpc = import_mpctools();

% Sizes and constants.
Delta = 0.1;
Nt = 50;

Nx = 3;
Nu = 2;
Ny = 2;
Nw = Nx;
Nv = Ny;

Ac = [-1, 1, 0; 0, -2, 2; 0, 0, -0.5];
Bc = [1, 0; 0, 1; 1, 1];
C = [1, 1, 0; 0, 1, 1];
Gc = Ac/(expm(Ac*Delta) - eye(Nx)); % Gives G = I in discrete time.

funcs = struct();

fc = @(x, u, w) Ac*x + Bc*u + Gc*w; % Continuous-time evolution.
funcs.fc = mpc.getCasadiFunc(fc, [Nx, Nu, Nw], {'x', 'u', 'w'}, {'fc'});

lin = mpc.getLinearizedModel(funcs.fc, ...
                             {zeros(Nx, 1), zeros(Nu, 1), zeros(Nw, 1)}, ...
                             {'A', 'B', 'G'}, 'Delta', Delta);

f = @(x, u, w) lin.A*x + lin.B*u + lin.G*w;
funcs.f = mpc.getCasadiFunc(f, [Nx, Nu, Nw], {'x', 'u', 'w'}, {'f'});

h = @(x) C*x;
funcs.h = mpc.getCasadiFunc(h, [Nx], {'x'}, {'h'});

% Pick noise covariances.
Q = 0.01*diag([0.1, 0.25, 0.05]);
[Qinv, Qhalf] = mpctools.spdinv(Q);

R = diag([0.5, 0.25]);
[Rinv, Rhalf] = mpctools.spdinv(R);

% Simulate noisy system.
rand('seed', 0);
x0 = [1; 2; 3];
omega = 2*pi()/(Nt*Delta);
t = (0:Nt)*Delta;
u = [sin(omega*t(1:end-1)); cos(omega*t(1:end-1))];
w = Qhalf*randn(Nw, Nt);
v = Rhalf*randn(Nv, Nt + 1);

x = NaN(Nx, Nt + 1);
x(:,1) = x0;
y = NaN(Ny, Nt + 1);

for k = 1:(Nt + 1)
    % Take measurement.
    y(:,k) = h(x(:,k)) + v(:,k);
    
    % Advance x.
    if k <= Nt
        x(:,k + 1) = f(x(:,k), u(:,k), w(:,k));
    end
end
mpctools.mpcplot(y, u, 'xnames', @(n) sprintf('y_{%d}', n));

% Define MHE stage costs and other arguments.
l = @(w, v) w'*Qinv*w + v'*Rinv*v;
funcs.l = mpc.getCasadiFunc(l, [Nw, Nv], {'w', 'v'}, {'l'});

lx = @(x, x0bar) 100*(x - x0bar)'*(x - x0bar);
funcs.lx = mpc.getCasadiFunc(lx, [Nx, Nx], {'x', 'x0bar'}, {'lx'});

N = struct('x', Nx, 'u', Nu, 'y', Ny, 't', Nt);

x0bar = [1.1; 1.9; 3.5];

solver = mpc.nmhe(funcs.f, funcs.l, funcs.h, u, y, N, funcs.lx, x0bar, ...
                  'Delta', Delta, 'verbosity', 5);

% Solve optimization and make a plot.
solver.solve();
xhat = solver.var.x;
vhat = solver.var.v;
err = C*xhat + vhat - y;
figure();
for i = 1:Nx
    subplot(Nx, 1, i);
    plot(t, xhat(i,:), '-xr', t, x(i,:), '-og', 0, x0bar(i), 'sk');
    ylabel(sprintf('x_{%d}', i));
    legend('Estimated', 'Actual', 'Prior');
end

