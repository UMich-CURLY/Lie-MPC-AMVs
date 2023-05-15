% Lotka-Volterra Fishing Problem
mpc = import_mpctools();

% System parameters.
Nt = 40;
Nx = 2;
Nu = 1;
Delta = 0.25;

function dxdt = ode(x, u)
    % ODE right-hand side.
    c1 = 0.4;
    c2 = 0.2;
    dxdt = [ ...
        x(1) - x(1)*x(2) - c1*x(1)*u(1);
        -x(2) + x(1)*x(2) - c2*x(2)*u(1);
    ];
end%function
f = mpc.getCasadiFunc(@ode, [Nx, Nu], {'x', 'u'}, {'f'}, 'rk4', true(), ...
                      'Delta', Delta, 'M', 4);

% Initial conditions, bounds, guess, etc.
x0 = [0.5; 0.7];
x = NaN(Nx, Nt + 1);
u = zeros(Nu, Nt); % Default guess is no fishing.
x(:,1) = x0;
for t = 1:Nt
    x(:,t + 1) = full(f(x(:,t), u(:,t)));
end
guess = struct('x', x, 'u', u);
lb = struct('x', zeros(Nx, 1), 'u', 0);
ub = struct('x', 2*ones(Nx, 1), 'u', 1);
udiscrete = true(Nu, 1);

% Stage cost.
function cost = stagecost(x, u)
    % Quadratic stage cost with steady-state x = [1; 1], u = [0].
    cost = (x(1) - 1).^2 + (x(2) - 1).^2 + 0.1*u(1);
end%function
l = mpc.getCasadiFunc(@stagecost, [Nx, Nu], {'x', 'u'}, {'l'});

% Create controller.
N = struct('x', Nx, 'u', Nu, 't', Nt);
controller = mpc.nmpc('f', f, 'l', l, 'N', N, 'x0', x0, 'lb', lb, 'ub', ub, ...
                      'guess', guess, 'udiscrete', udiscrete, ...
                      'solver', 'bonmin');
controller.solve();

% Plot.
mpc.mpcplot(controller.var.x, controller.var.u, ...
            'xnames', {'Prey', 'Predator'}, 'unames', {'Fishing'});

