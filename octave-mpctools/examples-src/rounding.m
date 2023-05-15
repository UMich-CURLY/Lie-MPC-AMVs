% Example of cumulative rounding for a discrete-valued system.
mpc = import_mpctools();

% Choose system.
Nx = 2;
Nu = 1;

Ac = [0, 1; -0.5, -0.25];
Bc = [0; 1];
Delta = 0.25;

umin = 0;
umax = 1;

[A, B] = mpc.c2d(Delta, Ac, Bc, 'return', 'AB');


usp = 0.1;
xsp = (eye(Nx) - A)\B*usp;

% Make controller.
f = mpc.getCasadiFunc(@(x,u) A*x + B*u, [Nx, Nu], {'x', 'u'}, {'f'});
l = mpc.getCasadiFunc(@(x,u) sum((x - xsp).^2) + sum((u - usp).^2), ...
                      [Nx, Nu], {'x', 'u'}, {'l'});
N = struct('x', Nx, 'u', Nu, 't', 25);
lb = struct('u', umin);
ub = struct('u', umax);
x0 = zeros(Nx, 1);
kwargs = struct('f', f, 'l', l, 'N', N, 'lb', lb, 'ub', ub, 'x0', x0);
controller = mpc.nmpc('**', kwargs);

% Helper functions.
function [udiscr, offset] = roundu(ucont, offset, quantum)
    % Rounds the given u trajectory using cumulative rounding. u must be scaled
    % to be between 0 and 1, and quantum must be positive and less than one.
    if ~isscalar(quantum) || quantum <= 0 || quantum > 1
        error('Invalid value for quantum!');
    end
    ucont = ucont/quantum;
    Nmax = floor(1/quantum);
    udiscr = zeros(size(ucont));
    wantsofar = 0;
    getsofar = offset/quantum;
    for t = 1:length(ucont)
        wantsofar = wantsofar + ucont(t);
        udiscr(t) = min(round(wantsofar - getsofar), Nmax);
        getsofar = getsofar + udiscr(t);
    end
    offset = offset + (udiscr(1) - ucont(1))*quantum;
    udiscr = udiscr*quantum;
end%function

function x = sim(f, x0, u)
    % Simulates the system x^+ = f(x,u).
    Nt = length(u);
    x = zeros(length(x0), 1 + Nt);
    x(:,1) = x0;
    for t = 1:Nt
        x(:,t + 1) = f(x(:,t), u(:,t));
    end
end%function

Nsim = 150;

% Get a long open-loop trajectory.
kwargs.N.t = Nsim;
olcontroller = mpc.nmpc('**', kwargs);

f = @(x, u) A*x + B*u;
quantum = 0.5;

olcontroller.solve();
sol = olcontroller.var;
rsol = struct('t', sol.t);
rsol.u = roundu(sol.u, 0, quantum);
rsol.x = sim(f, sol.x(:,1), rsol.u);

% Get closed-loop trajectory with rounding.
x = NaN(Nx, Nsim + 1);
x(:,1) = x0;
u = NaN(Nu, Nsim);
quanterr = 0;
for t = 1:Nsim
    controller.fixvar('x', 1, x(:,t));
    controller.solve();
    if ~isequal(controller.status, 'Solve_Succeeded')
        warning('Solver failed at time %d!', t);
        break
    end
    [u(:,t), quanterr] = roundu(controller.var.u(:,1), quanterr, quantum);
    x(:,t + 1) = f(x(:,t), u(:,t));
end

% Make a plot.
kwargs = struct('fig', figure());
mpc.mpcplot(sol.x, sol.u, sol.t, 'color', 'r', 'legend', 'OL Continuous', ...
            '**', kwargs);
mpc.mpcplot(rsol.x, rsol.u, rsol.t, 'color', 'b', 'legend', 'OL Discrete', ...
            '**', kwargs);
mpc.mpcplot(x, u, sol.t, 'color', 'g', 'legend', 'CL Discrete', ...
            '**', kwargs);

