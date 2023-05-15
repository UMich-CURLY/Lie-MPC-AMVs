% VDP Oscillator using MPCTools.
mpc = import_mpctools();
Delta = 0.5;
Nt = 20; % Controller horizon.
Nx = 2;
Nu = 1;

%<<ENDCHUNK>>

% Define system model.
f = @(x, u) [(1 - x(2)^2)*x(1) ...
             - x(2) + u(1); x(1)];
l = @(x, u) x(1)^2 + x(2)^2 + u^2;

kwargs = struct('funcname', 'f', ...
                'rk4', true(), ...
                'Delta', Delta, ...
                'M', 4, ...
                'quad', l, ...
                'quadname', 'l');
[fcasadi, lcasadi] = ...
    mpc.getCasadiFunc(f, [Nx, Nu], ...
                      {'x', 'u'}, ...
                      '**', kwargs);

%<<ENDCHUNK>>

% Choose parameters.
ulb = -1;
uub = 1;
x0 = [0; 1];
xf = [0; 0];

%<<ENDCHUNK>>

% Formulate NLP.
N = struct('x', Nx, 'u', Nu, 't', Nt);
solver = mpc.nmpc('f', fcasadi, 'l', lcasadi, ...
                  'N', N, 'x0', x0, 'xf', xf, ...
                  'lb', struct('u', ulb), ...
                  'ub', struct('u', uub), ...
                  'verbosity', 5);

%<<ENDCHUNK>>

% Solve the OCP.
solver.solve();

%<<ENDCHUNK>>

% Plot the solution.
mpc.mpcplot(solver.var.x, solver.var.u, ...
            Delta*(0:Nt));

