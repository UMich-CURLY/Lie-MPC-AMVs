% Example of gain scheduling for tracking a periodic solution.
mpc = import_mpctools();

% Define system.
Ac = [0, 0, 1, 0;
      0, 0, 0, 1;
      -1, 0, 0, 0;
      0, -1, 0, 0];
Bc = [0, 0; 0, 0; 1, 0; 0, 1];
[Nx, Nu] = size(Bc);
Delta = pi()/32;

[A, B] = mpc.c2d(Delta, Ac, Bc, 'return', 'AB');
f = mpc.getCasadiFunc(@(x, u) A*x + B*u, [Nx, Nu], {'x', 'u'}, {'f'});

function y = measurement(x)
    % Measure distance from the unit circle in position and velocity space.
    y = [sqrt(x(1)^2 + x(2)^2) - 1; sqrt(x(3)^2 + x(4)^2) - 1];
end%function
h = mpc.getCasadiFunc(@(x) measurement(x), [Nx], {'x'}, {'f'});

% Generate a circular trajectory and linearize the system along it.
Nbar = 64;
xbar = NaN(Nx, Nbar + 1);
xbar(:,1) = [1; 0; 0; 1];
ubar = zeros(Nu, Nbar);
for i = 1:Nbar
    xbar(:,i + 1) = full(f(xbar(:,i), ubar(:,i)));
end

% Get a series of LQRs around the trajectory.
R = eye(Nu);
Kcirc = NaN(Nu, Nx, Nbar);
Qcirc = NaN(Nx, Nx, Nbar);
for i = 1:Nbar
    C = mpc.getLinearizedModel(h, {xbar(:,i)}, {'C'}, 'deal', true());
    Qcirc(:,:,i) = C'*C;
    Kcirc(:,:,i) = -dlqr(A, B, Qcirc(:,:,i), R);
end

% Simulate trajectory just using the given control laws.
Nt = 128;
x0 = zeros(Nx, 1);
xnom = NaN(Nx, Nt + 1);
xnom(:,1) = x0;
unom = NaN(Nu, Nt);
for i = 1:Nt
    ii = mod(i - 1, Nbar) + 1; % Index for reference trajectory.
    unom(:,i) = Kcirc(:,:,ii)*(xnom(:,i) - xbar(:,ii)) + ubar(:,ii);
    xnom(:,i + 1) = A*xnom(:,i) + B*unom(:,i);
end

% Now try with MPC.
usizes = {Nx, Nu, Nx, Nu, [Nu, Nx]};
uargs = {'x', 'u', 'xbar', 'ubar', 'K'};
function uact = localcontrol(x, u, xbar, ubar, K)
    % Applies local control law with u as a deviation from the plan.
    uact = u + ubar + K*(x - xbar);
end%function
ufunc = mpc.getCasadiFunc(@localcontrol, usizes, uargs, {'ufunc'});

% We're repurposing y as the actual control input, hence the odd model below.
faug = mpc.getCasadiFunc(@(x,y) A*x + B*y, [Nx, Nu], {'x', 'y'}, {'faug'});

function cost = stagecost(x, u, xbar, Q)
    % Penalize deviations from x setpoint and adjustment to nominal u.
    dx = x - xbar;
    cost = dx'*Q*dx + u'*u;
end%function
l = mpc.getCasadiFunc(@stagecost, {Nx, Nu, Nx, [Nx, Nx]}, ...
                      {'x', 'u', 'xbar', 'Q'}, {'l'});

function cost = termcost(x, xbar)
    % High penalty on deviation from setpoint.
    dx = x - xbar;
    cost = 100*dx'*dx;
end%function
Vf = mpc.getCasadiFunc(@termcost, [Nx, Nx], {'x', 'xbar'}, {'Vf'});

par = struct('xbar', xbar, 'ubar', ubar, 'K', Kcirc, 'Q', Qcirc);

lb = struct('y', -1*ones(Nu, Nt));
ub = struct('y', 1*ones(Nu, Nt));

N = struct('x', Nx, 'u', Nu, 'y', Nu, 't', Nt);

controller = mpc.nmpc('f', faug, 'l', l, 'h', ufunc, 'Vf', Vf, 'x0', x0, ...
                      'lb', lb, 'ub', ub, 'N', N, 'par', par, ...
                      'finaly', false());
controller.solve();
disp(controller.status);

x = controller.var.x;
u = controller.var.y;

tbar = mod(0:(Nt - 1), Nbar) + 1; % Need to repeat reference trajectory.
xbar = [xbar(:,tbar), xbar(:,1)];
ubar = ubar(:,tbar);

% Make a timeseries plot.
kwargs = struct('fig', figure());
mpc.mpcplot(x, u, 'color', 'black', 'legend', 'Adjusted', '**', kwargs);
mpc.mpcplot(xnom, unom, 'color', 'blue', 'linestyle', '--', ...
            'legend', 'Nominal', '**', kwargs);
mpc.mpcplot(xbar, ubar, 'color', 'red', 'legend', 'Setpoint', ...
            'linestyle', ':', '**', kwargs);

