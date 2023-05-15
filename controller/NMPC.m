function [u, d] = NMPC(x0, xd, xid, param, d)

if d.init == 0 % build a optimization problem and solve
    addpath('./casadi-linux-matlabR2014b-v3.5.5/'); 
    addpath('./octave-mpctools/');
    addpath('./GNC/');
    addpath('./HYDRO/');
    
    import casadi.*
    d = [];
    d.c.mpc = import_mpctools();
    d.p.x0 = zeros(12,1);
    x_otter = [x0(11:13); x0(8:10); x0(5:7) ;quat2eul(x0(1:4)', 'XYZ')'];
   
    d.p.x0 = x_otter;
    d = build_setup(d, param);
    eul_xyz = quat2eul(xd(:, 1:4), 'XYZ');
    x_ref = zeros(param.Nt, param.Nx);
    x_ref(:,1:3) = xid(:,4:6);
    x_ref(:,4:6) = xid(:,1:3);
    x_ref(:,7:9) = xd(:,5:7);
    x_ref(:,10:12) = eul_xyz;
    d = create_NMPC(d, x_ref);
    d = solve_NMPC(d,1);
    d.init = 1;

else % update parameters and solve
    x_otter = [x0(11:13); x0(8:10); x0(5:7) ;quat2eul(x0(1:4)', 'XYZ')'];
    d.p.x0 = x_otter;
    eul_xyz = quat2eul(xd(:, 1:4), 'XYZ');
    x_ref = zeros(param.Nt, param.Nx);
    x_ref(:,1:3) = xid(:,4:6);
    x_ref(:,4:6) = xid(:,1:3);
    x_ref(:,7:9) = xd(:,5:7);
    x_ref(:,10:12) = eul_xyz;
    d.c.solvers.NMPC.par.xsp = x_ref';
    d = solve_NMPC(d,1);

end

u = d.s.u(:,1);
end

function d = solve_NMPC(d,t)
    
    d.c.solvers.NMPC.fixvar('x', 1, d.p.x0(:,1));
    d.c.solvers.NMPC.solve();
    d.c.solvers.NMPC.saveguess();
    d.s.u(:,t) = d.c.solvers.NMPC.var.u(:,1);
    
end

function d = create_NMPC(d, x_ref)
    
    % import dynamics
    ode_casadi_NMPC = d.c.mpc.getCasadiFunc(...
        @otter, ...
        [d.p.n_x, d.p.n_u], ...
        {'x', 'u'});
    
    % discretize dynamics in time
    F = d.c.mpc.getCasadiFunc(...
        ode_casadi_NMPC, ...
        [d.p.n_x, d.p.n_u], ...
        {'x', 'u'}, ...
        'rk4', true(), ...
        'Delta', d.p.T);

    % define cost functions

    d.xsp = x_ref';
    l = @(x, u, xsp) (x(1:9)-xsp(1:9))'*d.p.Q(1:9,1:9)*(x(1:9)-xsp(1:9)) ...
        +atan2(sin(x(10)-xsp(10)), cos(x(10)-xsp(10)))^2*d.p.Q(10,10) ...
        +atan2(sin(x(11)-xsp(11)), cos(x(11)-xsp(11)))^2*d.p.Q(11,11) ...
        +atan2(sin(x(12)-xsp(12)), cos(x(12)-xsp(12)))^2*d.p.Q(12,12) ...
        + u'*d.p.R*u;

    Vf = @(x, xsp) (x(1:9)-xsp(1:9))'*d.p.P(1:9,1:9)*(x(1:9)-xsp(1:9)) ...
        +atan2(sin(x(10)-xsp(10)), cos(x(10)-xsp(10)))^2*d.p.P(10,10) ...
        +atan2(sin(x(11)-xsp(11)), cos(x(11)-xsp(11)))^2*d.p.P(11,11) ...
        +atan2(sin(x(12)-xsp(12)), cos(x(12)-xsp(12)))^2*d.p.P(12,12);

    lcasadi = d.c.mpc.getCasadiFunc(l, [d.p.n_x, d.p.n_u, d.p.n_x], {'x', 'u', 'xsp'}, {'l'});
    Vfcasadi = d.c.mpc.getCasadiFunc(Vf, [d.p.n_x, d.p.n_x], {'x', 'xsp'}, {'Vf'});

    % define NMPC arguments
    commonargs.l = lcasadi;
    commonargs.lb.x = d.p.x_min_v;
    commonargs.ub.x = d.p.x_max_v;
    commonargs.lb.u = d.p.u_min_v;
    commonargs.ub.u = d.p.u_max_v;
    commonargs.Vf = Vfcasadi;
    % commonargs.ef = ef;
    
    par = struct('xsp', d.xsp);

    % define NMPC problem dimensions
    N.x = d.p.n_x; % state dimension
    N.u = d.p.n_u; % control input dimension
    N.t = d.p.N_NMPC; % time dimension (i.e., prediction horizon)
    
    % create NMPC solver
    d.c.solvers.NMPC = d.c.mpc.nmpc(...
        'f', F, ... % dynamics (discrete-time)
        'N', N, ... % problem dimensions
        'Delta', d.p.T, ... % timestep
        'timelimit', d.p.cal_time, ... % solver time limit (in seconds)
        'par', par, ...
        '**', commonargs); % arguments
    
end


function d= build_setup(d, param)
    % state constraints
    d.p.x_min = -Inf;
    d.p.x_max = Inf;
    
    % control input constraints
    d.p.u_min = param.umin(1);
    d.p.u_max = param.umax(1);
    
    % NMPC prediction horizon (in number of time steps)
    d.p.N_NMPC = param.Nt;
    
    % sampling time (in time units)
    d.p.T = param.dt;
    
    % number of state variables
    d.p.n_x = param.Nx;
    
    % number of control inputs
    d.p.n_u = param.Nu;
    
    % simulation length (in number of time steps)
    d.p.t_final = 2;
    
    % pre-allocate memory
    d.s.x = NaN(d.p.n_x,d.p.t_final);
    d.s.u = NaN(d.p.n_u,d.p.t_final);
    
    % set initial state
    % (d.p.x0 is an argument of the main function)
    d.s.x(:,1) = d.p.x0;
    
    % state constraints vector
    d.p.x_min_v = d.p.x_min*ones(d.p.n_x,1);
    d.p.x_max_v = d.p.x_max*ones(d.p.n_x,1);
    
    % control input constraints vector
    d.p.u_min_v = d.p.u_min*ones(d.p.n_u,1);
    d.p.u_max_v = d.p.u_max*ones(d.p.n_u,1);
    
    % weighting matrices of stage cost
    % d.p.P = param.P;
    d.p.P = diag(zeros(1,12));
    d.p.P(7:9,7:9) = param.P(4:6,4:6);
    d.p.P(10:12,10:12) = param.P(1:3,1:3);
    d.p.P(4:6,4:6) = param.P(7:9,7:9);
    d.p.P(1:3,1:3) = param.P(10:12,10:12);
    % d.p.Q = param.Q;
    d.p.Q = diag(zeros(1,12));
    d.p.Q(7:9,7:9) = param.Q(4:6,4:6);
    d.p.Q(10:12,10:12) = param.Q(1:3,1:3);
    d.p.Q(4:6,4:6) = param.Q(7:9,7:9);
    d.p.Q(1:3,1:3) = param.Q(10:12,10:12);
    d.p.R = param.R;
    
    d.p.alpha = 0;
    d.p.cal_time = param.mpc_cal_t;

end
