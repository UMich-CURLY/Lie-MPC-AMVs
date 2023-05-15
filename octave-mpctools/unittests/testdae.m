% Tests a simple DAE to make sure it matches its exact ODE.
mpc = import_mpctools();

% Define model.
Delta = 0.5;
Nx = 2;
Nz = 1;
Nu = 1;

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

daes = struct('ode', mpc.getCasadiDAE(Delta, fode));
solvers = {'idas', 'collocation'};
for i = 1:length(solvers)
    s = solvers{i};
    daes.(s) = mpc.getCasadiDAE(Delta, fdae, gdae, 'solver', s);
end


% Simulate.
Nsim = 100;
x0 = [1; 0];
z0 = zfunc(x0);
t = Delta*(0:Nsim);
u = cos(2*pi()*t/20);
x = struct();
z = struct();
solvers = [{'ode'}, solvers];
for i = 1:length(solvers)
    s = solvers{i};
    x.(s) = [x0, NaN(Nx, Nsim)];
    z.(s) = [z0, NaN(Nz, Nsim)];
    dae = daes.(s);
    for k = 1:Nsim
        if isequal(s, 'ode')
            xp = dae(x.(s)(:,k), u(:,k));
            zp = zfunc(xp);
        else
            [xp, zp] = dae(x.(s)(:,k), z.(s)(:,k), u(:,k));
        end
        x.(s)(:,k + 1) = full(xp);
        z.(s)(:,k + 1) = full(zp);
    end
    if ~isequal(s, 'ode')
        err = [x.(s) - x.ode; z.(s) - z.ode];
        fprintf('%s max error: %g\n', s, max(abs(err(:))));
    end
end

% Make a plot.
plotoptions = struct('fig', figure(), 'umarker', 'yes', 'unames', {{'z'}});
colors = 'rgbcmyk';
markers = 'os^dv68';
for i = 1:length(solvers)
    s = solvers{i};
    mpc.mpcplot(x.(s), z.(s), t, 'color', colors(i), 'marker', markers(i), ...
                'legend', s, '**', plotoptions);
end

% Check how close.
assert(x.ode, x.idas, 1e-4);
assert(z.ode, z.idas, 1e-4);
assert(x.ode, x.collocation, 1e-3); % Collocation is less accurate.
assert(z.ode, z.collocation, 1e-3);

