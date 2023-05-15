% Example rate-of-change penalty and constraints.
mpc = import_mpctools();

% Model.
A = [0.5, 1; 0, 2];
B = [0; 1];

Nx = size(A, 1);
Nu = size(B, 2);
Nt = 10;
Nsim = 20;

N = struct('x', Nx, 'u', Nu, 't', Nt);

f = mpc.getCasadiFunc(@(x,u) A*x + B*u, [Nx, Nu], {'x', 'u'}, {'f'});

% Stage cost and terminal penalty.
l = mpc.getCasadiFunc(@(x,Du) x'*x + 100*Du'*Du, [Nx, Nu], {'x', 'Du'}, {'l'});
Vf = mpc.getCasadiFunc(@(x) 100*x'*x, [Nx], {'x'}, {'Vf'});

% Bounds.
umax = 5;
Dumax = 1;

lb = struct();
lb.u = -umax;
lb.Du = -Dumax;

ub = struct();
ub.u = umax;
ub.Du = Dumax;

% Build controller.
tic();
controller = mpc.nmpc('f', f, 'l', l, 'Vf', Vf, 'N', N, 'lb', lb, 'ub', ub, ...
                    'uprev', zeros(Nu, 1));
fprintf('Building controller took %.5g s.\n', toc());

% Simulate.
x = NaN(Nx, Nsim + 1);
x(:,1) = ones(Nx, 1);
u = NaN(Nu, Nsim);
for t = 1:Nsim
    controller.fixvar('x', 1, x(:,t));
    controller.solve();
    fprintf('Step %d: %s\n', t, controller.status);
    
    u(:,t) = controller.var.u(:,1);
    x(:,t + 1) = controller.var.x(:,2);
    
    controller.par.uprev = u(:,t);
    controller.saveguess();
end

% Plot.
mpc.mpcplot(x, u);

