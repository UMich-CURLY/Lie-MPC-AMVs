% Tests inferring function type from a Casadi object.
asserttype = @(f, t) assert(mpctools.getfunctiontype_(f), t);

f = mpctools.getCasadiFunc(@(x) x.^2, 2, {'x'});
asserttype(f, 'SX');

% Note: mpower(x, 2) doesn't work for Casadi < 3.3.0.
f = mpctools.getCasadiFunc(@(x) x*x, {[2, 2]}, {'x'}, 'casaditype', 'MX');
asserttype(f, 'MX');

Delta = 1;
f = mpctools.getCasadiIntegrator(@(x) -x, Delta, 2, {'x'}, 'casaditype', 'SX');
asserttype(f, 'MX'); % Integrators are always MX, even when the underlying
                     % expression is SX.

f = mpctools.getCasadiIntegrator(@(x) -x, Delta, 2, {'x'}, 'casaditype', 'MX');
asserttype(f, 'MX');

