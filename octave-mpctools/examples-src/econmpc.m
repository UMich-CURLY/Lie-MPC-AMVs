% Example from "A Lyapunov Function for Economic Optimizing Model Predictive
% Control" by Diehl, Amrit, and Rawlings (IEEE Trans. Auto. Cont., 56(3)), 2011
mpc = import_mpctools();

% Sizes.
Nx = 2;
Nu = 1;
Nt = 30;
Nsim = 9;
Nc = 4;
Nrho = Nx + Nu;

% Parameters.
pars = struct();
pars.cAf = 1;
pars.cBf = 0;
pars.Vr = 10;
pars.kr = 1.2;
pars.cAs = 0.5;
pars.cBs = 0.5;
pars.Qs = 12;
Delta = 0.5;
Qmax = 20;
data = struct('Delta', Delta);

xs = [pars.cAs; pars.cBs];
us = [pars.Qs];

function dxdt = cstrmodel(x, u, pars)
    % Nonlinear CSTR model.
    Q = u(1);
    cA = x(1);
    cB = x(2);
    dxdt = [ ...
        Q/pars.Vr*(pars.cAf - cA) - pars.kr*cA;
        Q/pars.Vr*(pars.cBf - cB) + pars.kr*cA;
    ];
end%function
f = mpc.getCasadiFunc(@(x,u) cstrmodel(x, u, pars), [Nx, Nu], {'x', 'u'}, {'f'});

function cost = stagecost(x, u, rho, pars)
    % Economic cost with regularization.
    Q = u(1);
    cA = x(1);
    cB = x(2);
    cost = -2*Q*cB + 0.5*Q + rho(1)*(cA - pars.cAs).^2 ...
           + rho(2)*(cB - pars.cBs).^2 + rho(3)*(Q - pars.Qs).^2;
end%function
l = mpc.getCasadiFunc(@(x,u,rho) stagecost(x, u, rho, pars), [Nx, Nu, Nrho], ...
                      {'x', 'u', 'rho'}, {'l'});

function cost = termcost(x)
    % Quadratic terminal cost.
    cost = 1e5*x'*x;
end%function
Vf = mpc.getCasadiFunc(@(x) termcost(x), [Nx], {'x'}, {'Vf'});

% Make a poor guess.
guess = struct();
guess.x = repmat(xs, 1, Nt + 1);
guess.xc = repmat(xs, 1, Nc, Nt);
guess.u = repmat(us, 1, Nt);

% Build controller.
kwargs = struct();
kwargs.f = f;
kwargs.l = l;
kwargs.Vf = Vf;
kwargs.N = struct('x', Nx, 'u', Nu, 'c', Nc, 't', Nt);
kwargs.Delta = Delta;
kwargs.lb = struct('x', zeros(Nx, Nt + 1), 'u', zeros(Nu, Nt));
kwargs.ub = struct('u', Qmax*ones(Nu, Nt));
kwargs.par = struct('rho', zeros(Nrho, 1));
kwargs.discretel = false();

tic();
controller = mpc.nmpc('**', kwargs);
fprintf('Building controller took %.4g s.\n', toc());

% Simulate with regularization and without.
tic();
t = linspace(0, 2*pi(), 10);
t = t(1:end-1);
x0 = [(1 + 0.7*cos(t))*pars.cAs; (1 + 0.7*sin(t))*pars.cBs];
colors = jet(length(t));
lam = [-10; -20];
rhos = {zeros(Nrho, 1), 0.505*ones(Nrho, 1)};
titles = {'Economic', 'Regularized'};
for i = 1:length(rhos)
    fprintf('** %s Objective\n', titles{i});
    rho = rhos{i};
    controller.par.rho = rho;
    figure();
    hold('on');
    trajectories = cell(size(x0, 2), 1);
    for j = 1:size(x0, 2)
        fprintf('  Initial Condition %d', j);
        controller.saveguess(guess); % Revert to default guess.
        x = NaN(Nx, Nsim + 1);
        x(:,1) = x0(:,j);
        xc = NaN(Nx, Nc, Nsim);
        V = NaN(Nsim);
        Vrot = NaN(Nsim);
        for t = 1:Nsim
            controller.fixvar('x', 1, x(:,t));
            controller.solve();
            if ~strcmp(controller.status, 'Solve_Succeeded')
                error('Solver failure for rho %d, x0 %d, t %d!', i, j, t);
            end
            fprintf('.')
            V(t) = full(controller.sol.f - Nt*stagecost(xs, us, rho, pars) ...
                        - termcost(controller.var.x(:,end)));
            Vrot(t) = V(t) - lam'*(x(:,t) - xs);
            x(:,t + 1) = controller.var.x(:,2);
            xc(:,:,t) = controller.var.xc(:,:,1);
            controller.saveguess();
        end
        fprintf('\n');
        X = mpc.smushcolloc(x, xc);
        plot(x(1,:), x(2,:), 'o', 'color', colors(j,:));
        plot(X(1,:), X(2,:), '-', 'color', colors(j,:));
        trajectories{j} = struct('x', x, 'X', X, 'V', V, 'Vrot', Vrot);
    end
    title(titles{i});
    xlabel('c_A');
    ylabel('c_B', 'rotation', 0);
    
    data.(titles{i}) = trajectories;
end
fprintf('Simulation and plotting took %.4g s.\n', toc());

% <--
% Save data.
save('-v7', 'econmpc.mat', '-struct', 'data');
% -->

