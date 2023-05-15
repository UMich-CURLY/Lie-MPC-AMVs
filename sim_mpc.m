function logger = sim_mpc(q0, p0, w0, v0, dt, Nsim, mode, param)
X_ref = param.X_ref;
xi_ref = param.xi_ref;

% design cost function parameters
param.P = diag([1e-5,1e-5, 0.5...
    10,10,1e-5,...
    1e-5,1e-5,1e-5,...
    1,1e-5,1e-5]) * 10;
param.Q = diag([1e-5,1e-5, 0.5,...
    10,10,1e-5,...
    1e-5,1e-5,1e-5,...
    1,1e-5,1e-5]) * 1;
param.R = eye(2) * 1e-5;

%%
% log the sates
x0 = [q0;p0;w0;v0];
X = x0;
U = [];
Err = [];
Err2 = [];

logger = struct();
solver = struct();
solver.init = 0;
logger.taken_time = 0;
for i = 1:Nsim - param.Nt
    X_ref_rt = X_ref(i,:);
    xi_ref_rt = xi_ref(i:i+param.Nt-1, :);
    tic;
    if mode == 1
        [u, solver] = eMPC(x0, X_ref_rt,  xi_ref_rt, param, solver);
    elseif mode == 2
        X_ref_rt = X_ref(i:i+param.Nt-1, :);
        [u, solver] = NMPC(x0, X_ref_rt, xi_ref_rt, param, solver);
        X_ref_rt = X_ref(i,:);
    elseif mode == 3
        X_ref_rt = X_ref(i:i+param.Nt-1, :);
        [u, solver] = NMPC_simple(x0, X_ref_rt, xi_ref_rt, param, solver);
        X_ref_rt = X_ref(i,:);
    else
        assert(1)
    end
    temp_taken_time = toc;
    logger.taken_time = logger.taken_time + temp_taken_time;
    
    x_otter = [x0(11:13); x0(8:10); x0(5:7) ;quat2eul(x0(1:4)', 'XYZ')'];

    % simulation frequency 80 Hz (4 times of control frequency)
    for j = 1:4
        x_otter = x_otter + otter_true(x_otter,u, param.V_c, 1)* param.dt/4;
    end

    x0 = [ eul2quat(x_otter(10:12)', 'XYZ')'; x_otter(7:9); x_otter(4:6); x_otter(1:3)];

    Xd = [quat2rotm(X_ref_rt(1:4)), X_ref_rt(5:7)';...
        0,0,0,1];
    X0 = [quat2rotm(x0(1:4)'), x0(5:7);...
        0,0,0,1];
    Xerr = logm(X0^-1 * Xd);
    XXerr = X0^-1 * Xd;
    x00 = x0(2:7);
    x00(1:3) = [Xerr(3,2); Xerr(1,3); Xerr(2,1)]; %% R error
    x00(4:6) = Xerr(1:3,4); %% p error

    X = [X, x0];
    U = [U, u];
    Err = [Err, x00];
    Err2 = [Err2; [sqrt(sum(sum(logm(XXerr(1:3,1:3)).^2))/2), XXerr(1:3,4)'] ];
end
logger.X = X;
logger.U = U;
logger.Err = Err;
logger.Err2 = Err2;
end