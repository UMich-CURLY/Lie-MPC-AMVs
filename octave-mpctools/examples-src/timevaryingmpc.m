% Example MPC with time-varying model stage cost.
mpc = import_mpctools();

% Define parameters and system model.
T = 1;
Delta = 0.01;
Nx = 1;
Nu = 1;
Nt = 250;

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
function cost = l(x, u, xsp, usp, Q, R)
    dx = x - xsp;
    du = u - usp;
    cost = dx'*Q*dx + du'*R*du;
endfunction
lcasadi = mpc.getCasadiFunc(@l, {Nx, Nu, Nx, Nu, [Nx, Nx], [Nu, Nu]}, ...
                            {'x', 'u', 'xsp', 'usp', 'Q', 'R'}, {'l'});

% Build bounds, parameters, and N.
lb = struct('u', -ones(Nu, Nt));
ub = struct('u', ones(Nu, Nt));
par = struct('xsp', xsp, 'usp', usp, 'Q', Q, 'R', R);

N = struct('x', Nx, 'u', Nu, 't', Nt);

x0 = -2;

% Build solver and optimize.
solver = mpc.nmpc('f', fcasadi, 'l', lcasadi, 'N', N, 'lb', lb, 'ub', ub, ...
                  'x0', x0, 'par', par, 'verbosity', 3);
solver.solve();

% Make a plot.
mpc.mpcplot('x', solver.var.x, 'u', solver.var.u, ...
            'xsp', solver.par.xsp, 'usp', solver.par.usp);
subplot(1, 2, 1);
legend('x', 'x_{sp}', 'Location', 'SouthEast');


