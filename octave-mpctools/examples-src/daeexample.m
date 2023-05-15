% Comparison of ODE and DAE formulations for Van der Pol oscillator.
mpc = import_mpctools();

% Define model.
Delta = 0.5;
Nsim = 20;
Nx = 2;
Nz = 1;
Nu = 1;
Nc = 3;
Nt = 25;

function z = zfunc(x)
    z = (1 - x(2).^2)*x(1);
end%function

function dxdt = odefunc(x, z, u)
    if isempty(z)
        z = zfunc(x);
    end
    dxdt = [z(1) - x(2) + u(1); x(1)];
end%function
ode = @(x, u) odefunc(x, [], u);
dae = @(x, z, u) odefunc(x, z, u);

vdp = mpc.getCasadiIntegrator(ode, Delta, [Nx, Nu], {'x', 'u'}, {'vdp'});

fode = mpc.getCasadiFunc(ode, [Nx, Nu], {'x', 'u'}, {'fode'});
fdae = mpc.getCasadiFunc(dae, [Nx, Nz, Nu], {'x', 'z', 'u'}, {'fdae'});
gdae = mpc.getCasadiFunc(@(x, z) z - zfunc(x), [Nx, Nz], {'x', 'z'}, {'gdae'});

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
N = struct('x', Nx, 'u', Nu, 'c', Nc, 't', Nt);
kwargs = struct('l', l, 'Vf', Vf, 'lb', lb, 'ub', ub, 'Delta', Delta, ...
                'x0', [1; 0], 'verbosity', 3);

solvers = struct();
solvers.ode = mpc.nmpc('f', fode, 'N', N, '**', kwargs);
N.z = Nz;
solvers.dae = mpc.nmpc('f', fdae, 'g', gdae, 'N', N, '**', kwargs);
solvernames = fieldnames(solvers);
colors = {'r', 'b'};
markers = {'o', 's'};
cmarkers = {'x', '+'};

% Solve both and plot.
fig = figure();
for i = 1:length(solvernames)
    name = solvernames{i};
    solver = solvers.(name);
    solver.solve();
    fprintf('%s: %s\n', name, solver.status);
    
    mpc.mpcplot(solver.var.x, solver.var.u, 'xc', solver.var.xc, 'fig', fig, ...
                'legend', name, 'color', colors{i}, 'marker', markers{i}, ...
                'umarker', markers{i}, 'collocmarker', cmarkers{i});
end

