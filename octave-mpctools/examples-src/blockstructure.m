% Builds an MPC problem and looks at the block structure.
mpc = import_mpctools();

% Nonlinear model.
Delta = 0.5;
Nsim = 20;
Nx = 2;
Nu = 1;
Nc = 3;
Nt = 10;

function dxdt = ode(x, u)
    dxdt = [(1 - x(2).^2)*x(1) - x(2) + u(1); x(1)];
endfunction
f = mpctools.getCasadiFunc(@ode, [Nx, Nu], {'x', 'u'}, {'f'});

function l = stagecost(x, u)
    l = x'*x + u'*u;
endfunction
l = mpctools.getCasadiFunc(@stagecost, [Nx, Nu], {'x', 'u'}, {'l'});

function Vf = termcost(x)
    Vf = 10*x'*x;
endfunction
Vf = mpctools.getCasadiFunc(@termcost, [Nx], {'x'}, {'Vf'});

% Set bounds and initial condition.
x0 = [1; 0];
lb = struct('u', -0.75*ones(Nu, Nt));
ub = struct('u', ones(Nu, Nt));

% Build solver and solve.
N = struct('x', Nx, 'u', Nu, 'c', Nc, 't', Nt);
solver = mpctools.nmpc('f', f, 'N', N, 'l', l, 'Vf', Vf, 'x0', x0, ...
                       'lb', lb, 'ub', ub, 'Delta', Delta);
solver.solve();
disp(solver.status);

% Build Lagrangian.
f = solver.nlp.f;
g = solver.nlp.g;
x = solver.nlp.x;
lamx = casadi.SX.sym('lamx', size(x));
lamg = casadi.SX.sym('lamg', size(g));
lagrangian = f + lamg'*g + lamx'*x;
var = [x; lamg; lamx];

hessian = lagrangian.hessian(var);
sparsity = hessian.sparsity();

H = zeros(size(hessian));
H(sparsity.find()) = 1;

% Plot the second derivative of the Lagrangian.
figure();
hold('on');
colormap('gray');
imagesc(~logical(H));
blocks = [0, size(x, 1), size(lamg, 1), size(lamx, 1)];
xplot = [0, size(H, 1)] + 0.5;
for i = 1:length(blocks)
    yplot = sum(blocks(1:i))*[1, 1] + 0.5;
    plot(xplot, yplot, '-k', yplot, xplot, '-k');
end
ticks = cumsum(blocks(1:(end - 1))) + 0.5*blocks(2:end);
ticklabels = {'x', 'lamg', 'lamx'};
set(gca(), 'xtick', ticks, 'xticklabel', ticklabels, 'xaxislocation', 'top', ...
    'ytick', ticks, 'yticklabel', ticklabels, 'ydir', 'reverse');
axis('equal');

