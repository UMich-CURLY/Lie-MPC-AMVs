% Rolling ball example.
mpc = import_mpctools();

% Some options.
movingHorizon = false();
terminalConstraint = false();
terminalCost = true();

% Sizes and bounds.
Nx = 4;
Nu = 2;
Nt = 25;
Nsim = 100;

umax = 1;
xmax = 2;
cushion = 0.1;

% Model.
Ac = zeros(Nx, Nx);
Ac(1, 3) = 1;
Ac(2, 4) = 1;

Bc = zeros(Nx, Nu);
Bc(3, 1) = 1;
Bc(4, 2) = 1;

% Discretize and make a function.
Delta = 0.05;
mats = mpctools.c2d(Delta, Ac, Bc);
A = mats.A;
B = mats.B;

f = mpctools.getCasadiFunc(@(x,u) A*x + B*u, [Nx, Nu], {'x', 'u'}, {'f'});

% Stage cost.
function cost = stagecost(x, u)
    % Quadratic stage cost.
    cost = x(1:2)'*x(1:2) + u'*u;
end%function
l = mpctools.getCasadiFunc(@stagecost, [Nx, Nu], {'x', 'u'}, {'l'});

function cost = termcost(x)
    % Quadratic terminal cost.
    cost = 1000*x(1:2)'*x(1:2);
end%function
Vf = mpctools.getCasadiFunc(@termcost, [Nx], {'x'}, {'Vf'});

% Terminal constraint.
rmin = 0.001;
function slack = termcon(x, rmin)
    % Terminal constraint to be close to the origin.
    slack = x(1).^2 + x(2).^2 - rmin.^2;
end%function
ef = mpctools.getCasadiFunc(@(x) termcon(x, rmin), [Nx], {'x'}, {'ef'});

% Define "holes".
r = 0.25;
m = 3;
centers = linspace(0, xmax, m + 1);
centers = 0.5*(centers(1:end-1) + centers(2:end));

circles = zeros(length(centers).^2, 3);
k = 0;
for i = 1:length(centers)
    for j = 1:length(centers)
        k = k + 1;
        circles(k, :) = [centers(i), centers(j), r];
    end
end

function slack = nlcon(x, c)
    % Nonlinear constraints.
    x1 = x(1);
    x2 = x(2);
    Nc = size(c, 1);
    slack = cell(Nc, 1);
    for i = 1:Nc
        c1 = c(i, 1);
        c2 = c(i, 2);
        r = c(i, 3);
        slack{i} = r.^2 - (x1 - c1).^2 - (x2 - c2).^2;
    end
    slack = vertcat(slack{:});
end%function
e = mpc.getCasadiFunc(@(x) nlcon(x, circles), [Nx], {'x'}, {'e'});

% Build optimizer.
if ~movingHorizon
    Nt = Nsim;
    Nsim = 1;
end

x0 = [xmax; xmax; 0; 0];

lb = struct();
lb.x = [-cushion; -cushion; -inf(); -inf()];
lb.u = -umax*ones(Nu, 1);

ub = struct();
ub.x = [xmax + cushion; xmax + cushion; inf(); inf()];
ub.u = umax*ones(Nu, 1);

guess = struct();
guess.x = repmat(x0, 1, Nt + 1);

N = struct('x', Nx, 'u', Nu, 't', Nt);

if movingHorizon
    verb = 0;
else
    verb = 5;
end

kwargs = struct('f', f, 'l', l, 'e', e, 'N', N, 'lb', lb, 'ub', ub, ...
                'guess', guess, 'maxiter', 5000,...
                'verbosity', verb, 'casaditype', 'MX');
if terminalConstraint
    kwargs.ef = ef;
end
if terminalCost
    kwargs.Vf = Vf;
end
tic();
controller = mpc.nmpc('**', kwargs);
fprintf('Building controller took %g s.\n', toc());

% Solve optimization.
x = NaN(Nx, Nsim + 1);
x(:,1) = x0;
u = NaN(Nu, Nsim);
tic();
for t = 1:Nsim
    controller.fixvar('x', 1, x(:,t));
    controller.solve();
    fprintf('%3d: %20s\n', t, controller.status);
    if ~strcmp(controller.status, 'Solve_Succeeded')
        warning('Failure at time %d!', t);
        break
    end
    
    if movingHorizon
        controller.saveguess();
        x(:,t + 1) = controller.var.x(:,2);
        u(:,t) = controller.var.u(:,1);
    else
        x = controller.var.x;
        u = controller.var.u;
    end
end
fprintf('Simulation took %g s.\n', toc());

% Make a plot.
figure();
hold('on');
for i = 1:size(circles, 1)
    th = circles(i,3)*exp(1i*linspace(0, 2*pi(), 255));
    fill(circles(i,1) + real(th), circles(i,2) + imag(th), 'r');
end
plot(x(1,:), x(2,:), '-k', x(1,1), x(2,1), 'ok', x(1,end), x(2,end), 'xk');
axis([0, xmax, 0, xmax] + 2*cushion*[-1, 1, -1, 1]);
xlabel('x coordinate');
ylabel('y coordinate');

% <--
% Save data.
data = struct('x', x, 'xmax', xmax, 'circles', circles);
save('-v7', 'ballmaze.mat', '-struct', 'data');
% -->

