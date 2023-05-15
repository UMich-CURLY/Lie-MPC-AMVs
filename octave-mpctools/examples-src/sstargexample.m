% Example of steady-state target finder to approximate sqrt(2).
mpc = import_mpctools();
N = struct('x', 2, 'u', 2, 'y', 1, 's', 2);

f = mpc.getCasadiFunc(@(x,u) x - u, [N.x, N.u], {'x', 'u'}, {'f'});
h = mpc.getCasadiFunc(@(x) x(1)^2/x(2)^2 - 2, [N.x], {'x'}, {'h'});
e = mpc.getCasadiFunc(@(y) [-y; y], [N.y], {'y'}, {'e'});
l = mpc.getCasadiFunc(@(x, s) log(x(2)) + 100*sum(s), [N.x, N.s], ...
                      {'x', 's'}, {'l'});

umax = 100;
umin = 1;
lb = struct('u', umin*ones(N.u, 1), 'x', ones(N.x, 1));
ub = struct('u', umax*ones(N.u, 1));

guess = struct('u', umin*ones(N.u, 1));
guess.x = guess.u;
guess.y = full(h(guess.x));
guess.s = [min(guess.y, 0); max(guess.y, 0)];

kwargs = struct('N', N, 'f', f, 'h', h, 'e', e, 'l', l, 'lb', lb, 'ub', ub, ...
                'guess', guess, 'discretef', false());

% First optimize with continuous variables.
sstarg = mpc.sstarg('**', kwargs);
sstarg.solve();
ipopt = sstarg.var;

% Now try discrete variables.
kwargs.udiscrete = true(N.u, 1);
kwargs.solver = 'bonmin';
sstarg = mpc.sstarg('**', kwargs);
sstarg.solve();
bonmin = sstarg.var;

% Use exhaustive search.
[num, den] = meshgrid(umin:umax, umin:umax);
err = (num./den).^2 - 2;
obj = log(den) + 100*abs(err);
[~, imin] = min(obj(:));
exact = struct('x', [num(imin); den(imin)], 's', err(imin));

% Print solutions.
printsol = @(sol) fprintf('%10g/%10g = %10g (error %10g)\n', sol.x(1), ...
                         sol.x(2), sol.x(1)/sol.x(2), sum(sol.s));

fprintf('Ipopt (continuous):\n');
printsol(ipopt);

fprintf('Bonmin (discrete):\n');
printsol(bonmin);

fprintf('Exact (discrete):\n');
printsol(exact);

