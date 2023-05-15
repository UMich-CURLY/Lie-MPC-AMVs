% Tests single- and multiple-shooting on a linear and nonlinear problem.
mpc = import_mpctools();

function [ss, ms] = docomparison(f, l, Vf, N, x0, Delta)
% Compares single shooting and multiple shooting.
narginchk(6, 6);
mpc = import_mpctools();

% Create Casadi functions.
if ~isnan(Delta)
    kwargs = struct('rk4', true(), 'Delta', Delta);
else
    kwargs = struct();
end
fcasadi = mpc.getCasadiFunc(f, [N.x, N.u], {'x', 'u'}, {'f'}, '**', kwargs);
lcasadi = mpc.getCasadiFunc(l, [N.x, N.u], {'x', 'u'}, {'l'});
Vfcasadi = mpc.getCasadiFunc(Vf, [N.x], {'x'}, {'Vf'});

% Get single and multiple shooting controllers.
kwargs = struct('f', fcasadi, 'N', N, 'x0', x0, 'l', lcasadi, 'Vf', Vfcasadi);
singleshooting = mpc.nmpc('singleshooting', true(), '**', kwargs);
multipleshooting = mpc.nmpc('**', kwargs);

% Solve OCPs.
singleshooting.solve();
ss = singleshooting.var;

multipleshooting.solve();
ms = multipleshooting.var;

% Compare via plot.
fig = figure();
t = (0:N.t);
if ~isnan(Delta)
    t = t*Delta;
end
kwargs = struct('t', t, 'fig', fig, 'umarker', 'yes');
mpc.mpcplot(ss.x, ss.u, 'color', 'red', 'marker', 'o', ...
            'legend', 'Single Shooting', '**', kwargs);
mpc.mpcplot(ms.x, ms.u, 'color', 'green', 'marker', 'x', ...
            'legend', 'Multiple Shooting', '**', kwargs);

end%function

% Linear example.
A = 1;
B = 1;
Q = 1;
R = 1;
x0 = 1;

pkg('load', 'control');
[K, P] = dlqr(A, B, Q, R);

f = @(x, u) A*x + B*u;
l = @(x, u) x'*Q*x + u'*Q*u;
Vf = @(x) x'*P*x;

N = struct('x', 1, 'u', 1, 't', 10);

[ss, ms] = docomparison(f, l, Vf, N, x0, NaN());
assert(ss.x, ms.x, 1e-4);
assert(ss.u, ms.u, 1e-4);

% Nonlinear example.
Delta = 0.5;
x0 = [0; 1];
N = struct('t', 20, 'x', 2, 'u', 1);

f = @(x, u) [(1 - x(2)^2)*x(1) - x(2) + u(1); x(1)];
l = @(x, u) x(1)^2 + x(2)^2 + u^2;
Vf = @(x) 100*x'*x;

[ss, ms] = docomparison(f, l, Vf, N, x0, Delta);
assert(ss.x, ms.x, 1e-4);
assert(ss.u, ms.u, 1e-4);

