% Example periodic optimization with demand charge.
mpc = import_mpctools();
rand('state', 0);

Nt = 24;
Nx = 1; % x = [storage tank level]
Nu = 2; % u = [direct purchase, storage withdrawal]
Nd = 1;

f = mpc.getCasadiFunc(@(x, u) x - u(2), [Nx, Nu], {'x', 'u'}, {'f'});
e = mpc.getCasadiFunc(@(u, umax, delta) [delta - u(1) - u(2); u(1) - umax], ...
                      [Nu, 1, 1], {'u', 'umax', 'delta'}, {'e'});
l = mpc.getCasadiFunc(@(u, rho) rho*u(1), [Nu, 1], {'u', 'rho'}, {'l'});
Vf = mpc.getCasadiFunc(@(umax, rhomax) rhomax*umax, [1, 1], ...
                       {'umax', 'rhomax'}, {'Vf'});

N = struct('x', Nx, 'u', Nu, 't', Nt);

lb = struct('u', [0; -inf()], 'x', 0);

t = 2*pi()*(0:(Nt - 1))/24;
par = struct('delta', 1 - cos(t) + rand(size(t)), ...
             'rho', 2 - cos(t) + rand(size(t)), 'rhomax', 0);
customvar = {'umax'};

controller = mpc.nmpc('f', f, 'l', l, 'Vf', Vf, 'e', e, 'N', N, ...
                      'lb', lb, 'par', par, 'customvar', customvar, ...
                      'customvar', customvar, 'periodic', true());

% Find optimal periodic solution for varying rhomax.
rhomaxes = [0, 0.1, 0.5, 1, 5, 10];
sols = cell(size(rhomaxes));
for i = 1:length(rhomaxes)
    controller.par.rhomax = rhomaxes(i);
    controller.solve();
    fprintf('rhomax = %g: %s\n', rhomaxes(i), controller.status);
    sols{i} = controller.var;
end

% Make a plot of parameters and optimal solutions.
t = 0:Nt;
figure();
subplot(2, 1, 1);
stairs(t, par.rho([1:end,end]), '-k');
ylabel('Price');
subplot(2, 1, 2);
stairs(t, par.delta([1:end,end]), '-k');
ylabel('Demand');
xlabel('Time');

style = struct('fig', figure(), 'xnames', {{'Storage'}}, ...
               'unames', {{'purchase', 'storage'}});
colors = jet(length(sols));
for i = 1:length(sols)
    mpc.mpcplot(sols{i}.x, sols{i}.u, 'color', colors(i,:), ...
                'legend', sprintf('\\rho_{max} = %g', rhomaxes(i)), ...
                '**', style);
end

