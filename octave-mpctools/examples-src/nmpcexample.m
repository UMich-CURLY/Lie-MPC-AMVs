% Linear and nonlinear control of a system.
mpc = import_mpctools();

% Constants.
Delta = 0.1;
Nx = 2;
Nu = 2;

function dxdt = ode(x, u)
    % Nonlinear ode.
    dxdt = [-x(1) - (1 + u(1))*x(2); (1 + u(1))*x(1) + x(2) + u(2)];
endfunction

% Simulator function.
model = mpc.getCasadiIntegrator(@ode, Delta, [Nx, Nu], {'x', 'u'});

% Ger nonlinear model and also linear approximation.
ode_casadi = mpc.getCasadiFunc(@ode, [Nx, Nu], {'x', 'u'});
xss = zeros(Nx, 1);
uss = zeros(Nu, 1);
lin = mpc.getLinearizedModel(ode_casadi, {xss, uss}, {'A', 'B'}, Delta);

% Define stage cost and terminal weight.
function l = stagecost(x, u)
    % Quadratic stage cost.
    l = 100*x'*x + u'*u;
endfunction
l = mpc.getCasadiFunc(@stagecost, [Nx, Nu], {'x', 'u'}, {'l'});

function Vf = termcost(x)
    % Quadratic terminal penalty.
    Vf = 1000*x'*x;
endfunction
Vf = mpc.getCasadiFunc(@termcost, [Nx], {'x'}, {'Vf'});

% Define linear model.
A = lin.A;
B = lin.B;
Ffunc = @(x, u) A*x + B*u;
F = mpc.getCasadiFunc(Ffunc, [Nx, Nu], {'x', 'u'}, {'Flin'});

% Make controllers.
Nt = 25;
x0 = [2; 2];
commonargs = struct('l', l, 'Vf', Vf, ...
                    'lb', struct('u', -ones(Nu, Nt)), ...
                    'ub', struct('u', ones(Nu, Nt)));
Nlin = struct('x', Nx, 'u', Nu, 't', Nt);
Nnonlin = Nlin;
Nnonlin.c = 2; % Use collocation for nonlinear MPC.

solvers = struct();
solvers.LMPC = mpc.nmpc('f', F, 'N', Nlin, '**', commonargs);
solvers.NMPC = mpc.nmpc('f', ode_casadi, 'N', Nnonlin, 'Delta', Delta, ...
                        '**', commonargs);

% Simulate.
Nsim = 100;
x = struct();
u = struct();
controllers = fieldnames(solvers);
for i = 1:length(controllers)
    c = controllers{i};
    x.(c) = NaN(Nx, Nsim + 1);
    x.(c)(:,1) = x0;
    u.(c) = NaN(Nu, Nsim);
    for t = 1:Nsim
        % Set initial condition and solve.
        solvers.(c).fixvar('x', 1, x.(c)(:,t));
        solvers.(c).solve();
        
        % Print status.
        fprintf('%5s %d: %s\n', c, t, solvers.(c).status);
        if ~isequal(solvers.(c).status, 'Solve_Succeeded')
            warning('%s failed at time %d!', c, t);
            break
        end
        solvers.(c).saveguess();
        u.(c)(:,t) = solvers.(c).var.u(:,1);
        x.(c)(:,t + 1) = full(model(x.(c)(:,t), u.(c)(:,t)));
        
        % Stop early if the system is near the origin.
        if norm(x.(c)(:,t + 1)) < 1e-2
            fprintf('%s at origin after %d iterations\n', c, t);
            x.(c)(:,(t + 2):end) = 0;
            u.(c)(:,(t + 1):end) = 0;
            break
        end
    end
end

% Make timeseries and phase plots.
figure();
colors = struct('LMPC', 'r', 'NMPC', 'g');
for i = 1:length(controllers)
    c = controllers{i};
    color = colors.(c);
    mpc.mpcplot(x.(c), u.(c), 'color', color, 'fig', gcf());
end

figure();
hold('on');
for i = 1:length(controllers)
    c = controllers{i};
    color = colors.(c);
    plot(x.(c)(1,:), x.(c)(2,:), 'color', color);
end
xlabel('x_1');
ylabel('x_2');
legend(controllers{:}, 'Location', 'SouthEast');

