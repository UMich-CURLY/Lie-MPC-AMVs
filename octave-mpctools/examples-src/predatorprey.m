% Estimation of predator-prey dynamics.
mpc = import_mpctools();
randn('state', 0);

% Global option.
global MEASURE_PREY
MEASURE_PREY = false();

% Sizes.
Nx = 3;
Nu = 1;
Ny = 1 + double(MEASURE_PREY);
Delta = 0.1;

% Pick coefficients.
pars = struct();
pars.a = 0.5; % Kill rate for prey.
pars.b = 1; % Birth rate for prey.
pars.c = 0.5; % Birth rate for predators.
pars.d = 1; % Death rate for predators.
pars.N0prey = 500; % Scale factor for prey.
pars.N0pred = 250; % Scale factor for predators.

function dxdt = ode_pars(x, u, pars)
    % Continuous-time dynamics.
    Nprey = x(1)/pars.N0prey;
    Ntag = x(2)/pars.N0prey;
    Npred = x(3)/pars.N0pred;

    Rtag = u(1);

    dxdt = [
        pars.N0prey*(-pars.a*Npred*Nprey + pars.b*Nprey - Rtag*Nprey);
        pars.N0prey*(-pars.a*Npred*Ntag + Rtag*Nprey);
        pars.N0pred*(pars.c*Npred*(Nprey + Ntag) - pars.d*Npred);
    ];
end%function
ode = @(x, u) ode_pars(x, u, pars);

function y = measurement(x)
    % Measure fraction of tagged animals.
    global MEASURE_PREY
    Nprey = x(1);
    Ntag = x(2);
    Ntot = Nprey + Ntag;
    tagfrac = Ntag/Ntot;
    if MEASURE_PREY
        y = [tagfrac; Ntot];
    else
        y = tagfrac;
    end 
end%function

% Convert to Casadi functions and simulator.    
model = mpc.getCasadiIntegrator(ode, Delta, [Nx,Nu], {'x','u'});
f = mpc.getCasadiFunc(ode, [Nx, Nu], {'x','u'}, {'f'}, ...
                      'rk4', true(), 'Delta', Delta);
h = mpc.getCasadiFunc(@measurement, [Nx], {'x'}, {'h'});
x0 = [pars.N0prey; 0; pars.N0pred];
x0bar = [1.5*pars.N0prey; 0; 0.75*pars.N0pred]; % Initial estimate.

% Simulate dynamics.
Nsim = 250;
x = NaN(Nx, Nsim + 1);
x(:,1) = x0;
u = zeros(Nu, Nsim);
t = Delta*(0:(Nsim - 1));
u(1,:) = 0.1*(1 + sin(2*pi()*t/10));
y = NaN(Ny, Nsim + 1);
yclean = NaN(Ny, Nsim + 1);

w = 1 + 0.05*randn(Nx, Nsim + 1);
v = 1 + 0.01*randn(Ny, Nsim + 1);
for t = 1:(Nsim + 1)
    % Round x and take measurement.
    x(:,t) = max(round(x(:,t).*w(:,t)), 0);
    yclean(:,t) = measurement(x(:,t));
    y(:,t) = yclean(:,t).*v(:,t);
    
    % Simulate step.
    if t <= Nsim
        x(:,t + 1) = full(model(x(:,t), u(:,t)));
    end
end

% Now try MHE.
function cost = stagecost(w, v)
    global MEASURE_PREY
    cost = 0.1*w'*w;
    if MEASURE_PREY
        cost = cost + 100*v(1)^2 + 0.25*v(2)^2;
    else
        cost = cost + 1000*v(1)^2;    
    end    
end%function

function cost = prior(x, x0bar)
    dx = x - x0bar;
    cost = 0.01*dx'*dx;
end%function

l = mpc.getCasadiFunc(@stagecost, [Nx,Ny], {'w','v'});
lx = mpc.getCasadiFunc(@prior, [Nx,Nx], {'x','x0bar'});

Nt = 35; % Window size for MHE.
N = struct('x', Nx, 'u', Nu, 'y', Ny, 'w', Nx, 'v', Ny, 't', Nt);
guess = struct('x', repmat(x0, 1, Nt + 1));
lb = struct('x', zeros(Nx, Nt + 1));
mhe = mpc.nmhe('f', f, 'h', h, 'u', u(:,1:Nt), 'y', y(:,1:(Nt + 1)), 'l', l, ...
               'lx', lx, 'N', N, 'x0bar', x0bar, 'lb', lb, 'guess', guess, ...
               'wadditive', true());

% Simulate closed loop.
xhat = NaN(Nx, Nsim + 1);
for t = 1:(Nsim + 1)
    % Solve current MHE problem.
    mhe.solve()
    if mod(t, 25) == 1 || ~isequal(mhe.status, 'Solve_Succeeded')
        fprintf('Step %d of %d: %s\n', t, Nsim - Nt, mhe.status);
    end
    
    % If at end, save the remaining trajectory. Otherwise, cycle.
    if t + Nt > Nsim
        xhat(:,t:end) = mhe.var.x;
        break
    else
        xhat(:,t) = mhe.var.x(:,1);
        mhe.newmeasurement(y(:,t + Nt), u(:,t + Nt), mhe.var.x(:,2));
        mhe.saveguess();
    end
end
yhat = NaN(Ny, Nsim + 1);
for t = 1:(Nsim + 1)
    yhat(:,t) = measurement(xhat(:,t));
end    

% Make a plot.
function doplot(t, x, xhat, y, yhat, yclean)
    % Plot actual and estimated 
    xlabels = {'Untagged Prey', 'Tagged Prey', 'Predators'};
    ylabels = {'Total Prey', 'Tag Fraction'};
    figure();
    Nx = size(x, 1);
    Ny = size(y, 1);
    Nplots = Nx + Ny;
    for i = 1:Nx
        subplot(Nplots, 1, i);
        hold('on');
        plot(t, x(i,:), '-g', 'DisplayName', 'Actual');
        plot(t, xhat(i,:), '-r', 'DisplayName', 'Estimated');
        ylabel(xlabels{i});
        if i == 1
            legend('location', 'NorthWest');
        end
    end
    for i = 1:Ny
        subplot(Nplots, 1, i + Nx)
        hold('on');
        plot(t, y(i,:), '-b', 'DisplayName', 'Measured');
        plot(t, yclean(i,:), '-c', 'DisplayName', 'Noise-Free');
        plot(t, yhat(i,:), '-m', 'DisplayName', 'Estimated');
        ylabel(ylabels{i});
        if i == 1
            legend('location', 'NorthWest');
        end
    end
end%function
t = (0:Nsim)*Delta;
doplot(t, x, xhat, y, yhat, yclean);
