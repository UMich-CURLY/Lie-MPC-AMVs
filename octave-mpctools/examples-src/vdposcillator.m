% Van der Pol oscillator example.
mpc = import_mpctools();

% Choose collocation or rk4 for nonlinear MPC.
usecollocation = true();

% Define model.
Delta = 0.5;
Nsim = 20;
Nx = 2;
Nu = 1;
Nt = 10;

function dxdt = ode(x, u)
    dxdt = [(1 - x(2).^2)*x(1) - x(2) + u(1); x(1)];
endfunction

vdp = mpc.getCasadiIntegrator(@ode, Delta, [Nx, Nu], {'x', 'u'}, {'vdp'});
if usecollocation
    kwargs = struct('funcname', 'vdpode');
    Nc = 3;
else
    kwargs = struct('funcname', 'vdprk4', 'rk4', true(), 'Delta', Delta);
    Nc = 0;
end
fnonlin = mpc.getCasadiFunc(@ode, [Nx, Nu], {'x', 'u'}, '**', kwargs);

linmodel = mpc.getLinearizedModel(fnonlin, {zeros(Nx, 1), zeros(Nu, 1)}, ...
                                  {'A', 'B'}, Delta);
Flin = mpc.getCasadiFunc(@(x, u) linmodel.A*x + linmodel.B*u, [Nx, Nu], ...
                         {'x', 'u'}, {'vdplin'});

function l = stagecost(x, u)
    l = x'*x + u'*u;
endfunction
l = mpc.getCasadiFunc(@stagecost, [Nx, Nu], {'x', 'u'}, {'l'});

function Vf = termcost(x)
    Vf = 10*x'*x;
endfunction
Vf = mpc.getCasadiFunc(@termcost, [Nx], {'x'}, {'Vf'});

% Set bounds.
lb = struct('u', -0.75*ones(Nu, Nt));
ub = struct('u', ones(Nu, Nt));

% Build solvers.
N = struct('x', Nx, 'u', Nu, 't', Nt);
kwargs = struct('l', l, 'Vf', Vf, 'lb', lb, 'ub', ub);

solvers = struct();
solvers.LMPC = mpc.nmpc('f', Flin, 'N', N, '**', kwargs);
N.c = Nc; % Collocation for nonlinear mpc.
solvers.NMPC = mpc.nmpc('f', fnonlin, 'N', N, 'Delta', Delta, '**', kwargs);

% Simulate closed-loop.
data = struct();
controllers = fieldnames(solvers);
for i = 1:length(controllers)
    c = controllers{i};
    solver = solvers.(c);
    fprintf('Simulating %s\n', c);
    x = NaN(Nx, Nsim + 1);
    xc = NaN(Nx, Nc, Nsim);
    x(:,1) = [1; 0]; % Initial condition.
    u = NaN(Nu, Nsim);
    for k = 1:Nsim
        % Solve MPC problem.
        solver.fixvar('x', 1, x(:,k));
        solver.solve();
        if ~isequal(solver.status, 'Solve_Succeeded')
            warning('Solver failed at time %d!', t);
            break
        end
        
        % Simulate.
        u(:,k) = solver.var.u(:,1);
        x(:,k + 1) = full(vdp(x(:,k), u(:,k)));
    end
    
    % Store data.
    data.(c) = struct('x', x, 'u', u, 't', Delta*(0:Nsim));
end

% Make a phase plot.
colors = {'r', 'g'};
figure();
hold('on');
for i = 1:length(controllers)
    c = controllers{i};
    plot(data.(c).x(1,:), data.(c).x(2,:), ['-o', colors{i}]);
end
xlabel('Velocity');
ylabel('Position');
legend(controllers{:});

% Also make a timeseries plot for each.
for i = 1:length(controllers)
    c = controllers{i};
    mpc.mpcplot('x', data.(c).x, 'u', data.(c).u, 't', data.(c).t, ...
                'title', c, 'xnames', {'Velocity', 'Position'}, ...
                'unames', {'Force'});
end

% <--
% Save data.
save('-v7', 'vdposcillator.mat', '-struct', 'data');
% -->

