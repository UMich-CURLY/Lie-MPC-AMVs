f = @(x) -x;
Delta = 0.1;
F = mpctools.getCasadiFunc(f, [1], {'x'}, {'F'}, 'rk4', true(), ...
                           'Delta', Delta, 'M', 1);
integrator = mpctools.getCasadiIntegrator(f, Delta, [1], {'x'}, {'integrator'});

x0 = 1;
check1 = full(F(x0));
check2 = full(integrator(x0));
checkexact = exp(-Delta)*x0;

assert(check1, checkexact, 1e-5);
assert(check2, checkexact, 1e-5);

