% Example MPC using 1-norm penalty on x and u.
mpc = import_mpctools();

% Bad things happen to IPOPT if you don't reformulate.
reformulateAbsoluteValues = true();

% Model and sizes.
A = [0.5, 1; 0, 2];
B = [0; 1];

Nx = size(A, 2);
Nu = size(B, 2);
Nt = 10;
Nsim = 25;

N = struct('x', Nx, 'u', Nu, 't', Nt);

f = mpc.getCasadiFunc(@(x,u) A*x + B*u, [Nx, Nu], {'x', 'u'}, {'f'});

% Stage cost and terminal cost.
if reformulateAbsoluteValues
    l = mpc.getCasadiFunc(@(absx, absu) sum(absx) + 10*sum(absu), [Nx, Nu], ...
                      {'absx', 'absu'}, {'l'});
    Vf = mpc.getCasadiFunc(@(absx) 100*sum(absx), [Nx], {'absx'}, {'Vf'});
else
    l = mpc.getCasadiFunc(@(x, u) sum(abs(x)) + 10*sum(abs(u)), [Nx, Nu], ...
                          {'x', 'u'}, {'l'});
    Vf = mpc.getCasadiFunc(@(x) 100*sum(abs(x)), [Nx], {'x'}, {'Vf'});
end

% Bounds.
lb = struct();
lb.u = -5;

ub = struct();
ub.u = 5;

% Build controller.
controller = mpc.nmpc('f', f, 'l', l, 'Vf', Vf, 'N', N, 'lb', lb, 'ub', ub);

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
    
    controller.saveguess();
end

% Plot.
mpc.mpcplot(x, u);

