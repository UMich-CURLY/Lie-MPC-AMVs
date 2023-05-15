% MHE on a spring pendulum using a DAE formulation.
mpc = import_mpctools();
randn('state', 0);

Nx = 4; % Bob position and velocity (2 dimensions)
Nz = 1; % Spring length.
Ny = 2; % Only measure position.

Nw = 2; % Disturbances on acceleration.
Nv = Ny;

% Specify gravity, mass, spring constant, and rest length.
pars = struct('g', -1, 'm', 1, 'k', 1, 'l0', 1);

function dxdt = ode(x, z, w, pars)
    T = pars.k*(z(1) - pars.l0); % Tension in spring.
    dxdt = [
        x(3);
        x(4);
        -T*x(1)/pars.m + w(1);
        pars.g - T*x(2)/pars.m + w(2);
    ];
end%function
f = mpc.getCasadiFunc(@(x, z, w) ode(x, z, w, pars), ...
                      [Nx, Nz, Nw], {'x', 'z', 'w'}, {'f'});

function err = algebra(x, z)
    err = x(1)^2 + x(2)^2 - z(1)^2;
end%function
g = mpc.getCasadiFunc(@(x, z) algebra(x, z), [Nx, Nz], {'x', 'z'}, {'g'});

function y = measurement(x)
    y = [x(1,:); x(2,:)]; % Only measure position.
end%function
h = mpc.getCasadiFunc(@(x) measurement(x), [Nx], {'x'}, {'h'});

% Simulate system.
Delta = 0.1;
dae = mpc.getCasadiDAE(Delta, f, g);

Nt = 100;
xsim = NaN(Nx, Nt + 1);
xsim(:,1) = [1; 1; 0; 0];
zsim = NaN(Nz, Nt + 1);
zsim(:,1) = hypot(xsim(1,1), xsim(2,1));

sig_w = 0.01;
sig_v = 0.001;

wsim = sig_w*randn(Nw, Nt);
vsim = sig_v*randn(Nv, Nt + 1);

for i = 1:Nt
    [xp, zp] = dae(xsim(:,i), zsim(:,i), wsim(:,i));
    xsim(:,i + 1) = full(xp);
    zsim(:,i + 1) = full(zp);
end
ysim = measurement(xsim) + vsim;

% Now try MHE to infer velocity and spring length using only position.
l = mpc.getCasadiFunc(@(w, v) w'*w/sig_w.^2 + v'*v/sig_v.^2, ...
                      [Nw, Nv], {'w', 'v'}, {'l'});

N = struct('x', Nx, 'z', Nz, 'y', Ny, 'w', Nw, 'v', Nv, 'c', 2, 't', Nt);

guess = struct();
guess.x = zeros(Nx, Nt + 1);
guess.x(1:2,:) = ysim;
guess.z = hypot(guess.x(1,:), guess.x(2,:));

fprintf('Building estimator...');
mhe = mpc.nmhe('f', f, 'g', g, 'h', h, 'l', l, 'y', ysim, 'N', N, ...
               'Delta', Delta, 'guess', guess, 'verbosity', 3);
mhe.solve();

kwargs = struct('fig', figure(), 't', (0:Nt)*Delta, ...
                'xnames', {{'x', 'y', 'v_x', 'v_y'}}, 'unames', {{'l'}});
mpc.mpcplot(xsim, zsim, 'marker', 'o', 'color', 'g', 'legend', 'actual', ...
            '**', kwargs);
mpc.mpcplot(mhe.var.x, mhe.var.z, 'marker', 's', 'color', 'r', ...
            'legend', 'estimated', '**', kwargs);

