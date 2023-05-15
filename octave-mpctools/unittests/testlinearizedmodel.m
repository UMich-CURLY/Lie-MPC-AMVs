% Tests linearization of a linear model.
A = [0.5, 1; 0, 0.5];
B = [0; 1];
F = mpctools.getCasadiFunc(@(x, u) A*x + B*u, [2, 1], {'x', 'u'}, {'F'});

x0 = [0; 0];
u0 = [0];
model = mpctools.getLinearizedModel(F, {x0, u0}, {'A', 'B'});

assert(A, model.A, 1e-10);
assert(B, model.B, 1e-10);

