% Simple system identification for a stratified tank using parest.
% **Note**: This script makes use of mpctools.parest, which is currently
% experimental. Thus, this example may change in the future along with parest.
mpc = import_mpctools();

pars = struct( ...
    'Vmax', 7.5, ... % Tank volume in m^3
    'rhocp', 4184, ... % Water volumetric heat capacity in kJ/m^3*K
    'Tchws', 5.5, ... % Chilled water supply temperature (*C)
    'Tchwr', 12.2, ... % Chilled water return temperature.
    'Delta', 1); % Timestep in hours.
pars.Tref = pars.Tchwr; % Reference temperature (*C)
kactual = 0.005*pars.Vmax; % Heat transfer coefficient between layers. (m^3/h)

function dxdt = tankmodel(x, u, k, pars)
    % Nonlinear dynamic model for storage tank. In terms of volumes and
    % enthalpies. For these purposes, H = rho*C*H*(T - T0) with T temperature.
    % u gives volumetric flow rates (enthalpies are calculated using pars).
    narginchk(4, 4);
    
    % Components of x.
    Hcold = x(1);
    Vcold = x(2);
    Hhot = x(3);
    Vhot = x(4);
    
    % Components of u.
    vplus = u(1);
    hplus = pars.rhocp*(pars.Tchws - pars.Tref);
    vminus = u(2);
    hminus = pars.rhocp*(pars.Tchwr - pars.Tref);
    
    % Nonlinear model. Need to be careful about division.
    hcold = (Hcold + hplus*eps())/(Vcold + eps());
    hhot = (Hhot + hminus*eps())/(Vhot + eps());
    dxdt = [
        hplus*vplus - hcold*vminus + k*(hhot - hcold);
        vplus - vminus;
        -hhot*vplus + hminus*vminus - k*(hhot - hcold);
        -vplus + vminus;
    ];
end%function

% Define functions.
N = struct('x', 4, 'u', 2, 'y', 4);
simulator = mpc.getCasadiIntegrator(@(x, u) tankmodel(x, u, kactual, pars), ...
                                    pars.Delta, [N.x, N.u], {'x', 'u'}, ...
                                    {'simulator'});
fnonlinear = mpc.getCasadiFunc(@(x, u, k) tankmodel(x, u, k, pars), ...
                               [N.x, N.u, length(kactual)], {'x', 'u', 'k'}, ...
                               {'fnonlinear'});
h = mpc.getCasadiFunc(@(x) x, [N.x], {'x'}, {'h'}); % Full state is measured.

% Choose initial condition and simulate.
Vcold0 = 0.01*pars.Vmax;
Hcold0 = pars.rhocp*Vcold0*(pars.Tchws - pars.Tref);
Vhot0 = pars.Vmax - Vcold0; 
Hhot0 = pars.rhocp*Vhot0*(pars.Tchwr - pars.Tref);

Nsim = 72;
xsim = NaN(N.x, Nsim + 1);
xsim(:,1) = [Hcold0; Vcold0; Hhot0; Vhot0];
u = 0.125*pars.Vmax*sin(2*pi()*(1:Nsim)*pars.Delta/24);
usim = [max(u, 0); -min(u, 0)]; % Split into charge and discharge flows.
ysim = NaN(N.y, Nsim + 1);
for t = 1:(Nsim + 1)
    ysim(:,t) = full(h(xsim(:,t)));
    if t <= Nsim
        xsim(:,t + 1) = full(simulator(xsim(:,t), usim(:,t)));
    end
end

% Build Nonlinear parameter estimator.
lb = struct('x', [-inf(); 0; -inf(); 0], 'k', 0);
guess = struct('x', xsim, 'k', 1);
N.t = Nsim;
N.c = 2; % Use collocation.
estimator = mpctools.parest('f', fnonlinear, 'h', h, 'par', {'k'}, ...
                            'data', struct('y', ysim, 'u', usim), 'N', N, ...
                            'lb', lb, 'guess', guess, 'Delta', pars.Delta, ...
                            'verbosity', 4);
estimator.solve()
xnonlinear = estimator.var.x;
knonlinear = estimator.var.k;

% Make a plot.
function tankplot(x, pars, lspecs, fig)
    % tankplot(x, pars, [lspecs={'-ob', '-sr'}], [fig])
    %
    % Plots the given tank trajectory. lspecs should give lspecs for plotting
    % the cold and hot states respectively. fig should be the figure handle to
    % use for the plots.
    narginchk(2, 4);
    if nargin() < 3
        lspecs = {'-ob', '-sr'};
    end
    if nargin() < 4
        fig = figure();
    end
    figure(fig);
    t = (0:(size(x, 2) - 1))*pars.Delta;
    hscale = 1e-3; % Convert from kWh to MWh.
    
    % Plot volume.
    subplot(3, 1, 1);
    hold('on');
    plot(t, x(2,:), lspecs{1}, t, x(4,:), lspecs{2});
    ylabel('Volume (m^3)');
    
    % Plot temperature.
    Tcold = pars.Tref + x(1,:)./(pars.rhocp*x(2,:));
    Thot = pars.Tref + x(3,:)./(pars.rhocp*x(4,:));
    subplot(3, 1, 2);
    hold('on');
    plot(t, Tcold, lspecs{1}, t, Thot, lspecs{2});
    ylabel('Temperature (*C)');
    
    % Plot enthalpy.
    subplot(3, 1, 3);
    hold('on');
    plot(t, hscale*x(1,:), lspecs{1}, t, hscale*x(3,:), lspecs{2});
    ylabel('Enthalpy (MWh)');
    xlabel('Time (h)');
end%function

fig = figure();
tankplot(xsim, pars, {'ob', 'sr'}, fig);
tankplot(xnonlinear, pars, {'-b', '-r'}, fig);

% Also try with simple linear model.
N = struct('x', 1, 'u', 1, 't', Nsim, 'y', 1);
data = struct('y', ysim(1,:), 'u', usim(1,:) - usim(2,:));
flinear = mpc.getCasadiFunc(@(x, u, a, b) a*x + b*u, [1, 1, 1, 1], ...
                            {'x', 'u', 'a', 'b'}, {'flinear'});
hlinear = mpc.getCasadiFunc(@(x) x, [1], {'x'}, {'hlinear'});
estimator = mpctools.parest('f', flinear, 'h', hlinear, 'par', {'a', 'b'}, ...
                            'data', data, 'N', N, 'verbosity', 4);
estimator.solve();
xlinear = [estimator.var.x; NaN(3, Nsim + 1)];
alinear = estimator.var.a;
blinear = estimator.var.b;

figure();
t = (0:Nsim)*pars.Delta;
scale = 1e-3;
plot(t, scale*xsim(1,:), 'ok', t, scale*xnonlinear(1,:), '-xg', ...
     t, scale*xlinear(1,:), '-+m');
legend('Data', 'Nonlinear Model', 'Linear Model', 'Location', 'SouthEast');
xlabel('Time (h)');
ylabel('Cold Enthalpy (MWh)');

