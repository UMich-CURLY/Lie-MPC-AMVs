% Test MPCTools on a simple linear MPC problem.
mpc = import_mpctools();
rand('state', 0);
randn('state', 0);

Nx = 2;
Nu = 1;
Ny = 1;
Nt = 5;
Nsim = 50;

% Choose system and inputs.
A = rand(Nx, Nx, Nsim);
B = rand(Nx, Nu, Nsim);
C = rand(Ny, Nx, Nsim + 1);

function H = spdrand(n, t, rho)
    % H = spdrand(n, t, [rho=0.1])
    %
    % Returns an array of t random n by n symmetric positive-definite matrices.
    % rho gives the multiple of the identity matrix to add to each matrix.
    narginchk(2, 3);
    if nargin() < 3
        rho = 0.1;
    end
    rho = rho*eye(n);
    H = NaN(n, n, t);
    for i = 1:t
        h = rand(n, n);
        H(:,:,i) = (h*h') + rho;
    end
end%function

Q = spdrand(Nx, Nsim);
R = spdrand(Ny, Nsim + 1);

u = sin(2*pi()*(0:(Nsim - 1))/25);
x = NaN(Nx, Nsim + 1);
x(:,1) = ones(Nx, 1);

w = randn(Nx, Nsim);
v = randn(Ny, Nsim + 1);

x0bar = x(:,1) + [50; -50];
P0 = eye(Nx);

% Simulate the system.
for i = 1:(Nsim + 1)
    y(:,i) = C(:,:,i)*x(:,i) + v(:,i);
    if i <= Nsim
        x(:,i + 1) = A(:,:,i)*x(:,i) + B(:,:,i)*u(:,i) + w(:,i);
    end
end

% Apply the time-varying Kalman Filter.
xkf = NaN(Nx, Nsim + 1); % This will be \hat{x}
xhatm = x0bar;
P = P0;
for i = 1:(Nsim + 1)
    L = (P*C(:,:,i)')/(C(:,:,i)*P*C(:,:,i)' + R(:,:,i));
    xkf(:,i) = xhatm + L*(y(:,i) - C(:,:,i)*xhatm); % Filter.
    
    P = P - L*C(:,:,i)*P;
    if i <= Nsim
        xhatm = A(:,:,i)*xkf(:,i) + B(:,:,i)*u(:,i);
        P = A(:,:,i)*P*A(:,:,i)' + Q(:,:,i);
    end
end

% Apply MHE.
f = mpc.getCasadiFunc(@(x, u, w, A, B) A*x + B*u + w, ...
                      {Nx, Nu, Nx, [Nx, Nx], [Nx, Nu]}, ...
                      {'x', 'u', 'w', 'A', 'B'}, {'f'});
h = mpc.getCasadiFunc(@(x, C) C*x, {Nx, [Ny, Nx]}, {'x', 'C'}, {'h'});
l = mpc.getCasadiFunc(@(w, v, Qinv, Rinv) w'*Qinv*w + v'*Rinv*v, ...
                      {Nx, Ny, [Nx, Nx], [Ny, Ny]}, ...
                      {'w', 'v', 'Qinv', 'Rinv'}, {'l'});

N = struct('x', Nx, 'u', Nu, 'y', Ny, 't', Nt);
par = struct('x0bar', x0bar, 'Qinv', mpc.spdinv(Q(:,:,1:Nt)), ...
             'Rinv', mpc.spdinv(R(:,:,1:(Nt + 1))), 'Pinv', mpc.spdinv(P0), ...
             'A', A(:,:,1:Nt), 'B', B(:,:,1:Nt), 'C', C(:,:,1:(Nt + 1)), ...
             'y', y(:,1:(Nt + 1)), 'u', u(:,1:Nt));
kwargs = struct('f', f, 'h', h, 'l', l, 'N', N, 'par', par, ...
                'isQP', true(), 'solver', 'ipopt');

% Simulate different prior updates.
priorupdates = {'smoothing', 'filtering'};
mhes = struct();
xmhes = struct();
for k = 1:length(priorupdates)
    key = priorupdates{k};
    fprintf('*** Simulating MHE with %s update.\n', key);
    mhe = mpc.nmhe('priorupdate', key, '**', kwargs);
    xmhe = NaN(Nx, Nsim + 1);
    for i = 0:Nsim
        if i <= Nt
            mhe.truncatehorizon(i);
        else
            mhe.cyclepar('Qinv', mpc.spdinv(Q(:,:,i)), ...
                         'Rinv', mpc.spdinv(R(:,:,i + 1)), ...
                         'A', A(:,:,i), 'B', B(:,:,i), 'C', C(:,:,i + 1), ...
                         'y', y(:,i + 1), 'u', u(:,i));
        end
        mhe.solve();
        fprintf('i = %d: %s\n', i, mhe.status);
        if ~isequal(mhe.status, 'Solve_Succeeded')
            warning('Solver failed!');
            break
        end
        mhe.saveestimate();
        xmhe(:,i + 1) = mhe.history(1).xhat(:,end);
    end
    
    mhes.(key) = mhe;
    xmhes.(key) = xmhe;
end

% Plot and compare.
fprintf('*** Plotting results.\n');
plotstyle = struct('fig', figure(), 'legendloc', 'NorthEast', ...
                   'unames', {repmat({'log_{10} |MHE - KF|'}, 1, 2)});
mpc.mpcplot(xkf, NaN(size(x)), 'legend', 'KF', 'color', 'k', 'marker', 'x', ...
            '**', plotstyle);
colors = {'m', 'c'};
markers = {'s', 'o'};
for k = 1:length(priorupdates)
    key = priorupdates{k};
    err = max(log10(abs(xkf - xmhes.(key))), -20);
    fprintf('%s max error: 10^%g\n', key, max(err(:)));
    mpc.mpcplot(xmhe, err, 'legend', sprintf('%s MHE', key), ...
                'color', colors{k}, 'marker', markers{k}, '**', plotstyle);
end
mpc.mpcplot(x, NaN(size(x)), 'legend', 'Actual', 'color', 'b', ...
            'linestyle', '--', '**', plotstyle);

