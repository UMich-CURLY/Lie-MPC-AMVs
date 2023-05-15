% Test MPCTools on a simple linear MPC problem.
mpc = import_mpctools();
randn('seed', 0);

% Define system.
A = 0.9;
B = 1;
G = 1;
C = 1;
Q = 0.01;
R = 0.1;
xhat0 = 1.5;
x0 = 2;

N = struct('x', 1, 'u', 1, 'y', 1, 'w', 1, 't', 20);

pkg('load','control');
[L, Pm] = dlqe(A, G, C, Q, R);
Pminv = mpctools.spdinv(Pm);

% Simulate system with noise.
[Qinv, Qhalf] = mpctools.spdinv(Q);
w = Qhalf*randn(N.w, N.t);

[Rinv, Rhalf] = mpctools.spdinv(R);
v = Rhalf*randn(N.y, N.t + 1);

u = 0.1*sin(2*pi()/N.t*(1:N.t));
x = NaN(N.x, N.t + 1);
y = NaN(N.y, N.t + 1);
x(:,1) = x0;
for i = 1:(N.t + 1)
    y(:,i) = C*x(:,i) + v(:,i);
    if i <= N.t
        x(:,i + 1) = A*x(:,i) + B*u(:,i) + G*w(:,i);
    end
end

% Define MHE problem.
f = @(x, u, w) A*x + B*u + G*w;
h = @(x) C*x;
l = @(w, v) w'*Qinv*w + v'*Rinv*v;
lx = @(x, x0bar) (x - x0bar)'*Pminv*(x - x0bar);

% Create Casadi functions.
fcasadi = mpc.getCasadiFunc(f, [N.x, N.u, N.w], {'x', 'u', 'w'}, {'f'});
hcasadi = mpc.getCasadiFunc(h, [N.x], {'x'}, {'h'});
lcasadi = mpc.getCasadiFunc(l, [N.w, N.y], {'w', 'v'}, {'l'});
lxcasadi = mpc.getCasadiFunc(lx, [N.x, N.x], {'x', 'x0bar'}, {'lx'});

solver = mpc.nmhe('f', fcasadi, 'N', N, 'x0bar', xhat0, 'l', lcasadi, ...
                  'lx', lxcasadi, 'h', hcasadi, 'u', u, 'y', y);

% Solve the full-information problem for t = 0 to Nt.
xmhe = NaN(N.x, N.t + 1, N.t + 1);
for t = 0:N.t
    solver.truncatehorizon(t);
    solver.solve();
    fprintf('t = %d: %s\n', t, solver.status);
    xmhe(:,:,t + 1) = solver.var.x;
end
plotinds = sub2ind(size(xmhe), ones(1, N.t + 1), 1:(N.t + 1), 1:(N.t + 1));

% Now use Kalman Filter.
xkf = NaN(N.x, N.t + 1);
xhat = xhat0;
for t = 1:(N.t + 1)
    xkf(:,t) = xhat + L*(y(:,t) - C*xhat);
    if t <= N.t
        xhat = A*xkf(:,t) + B*u(:,t);
    end
end

% Compare.
figure()
hold('on');
plot(0:N.t, xmhe(plotinds), '-om');
plot(0:N.t, xkf, '-xk');
plot(0:N.t, x, '-b');
legend('MHE', 'KF', 'Actual');
for t = 0:N.t
    plot(0:N.t, xmhe(:,:,t + 1), ':m');
end
ylabel('x', 'rotation', 0);
xlabel('Time');

