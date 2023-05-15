clear all; clc;close all;
% clearvars -except iii jjj kkk test_param1 test_param2 test_param3
addpath(genpath("./"));
%%
rng(2);
Ns = 10; 
dt = 0.05; % MPC simulation time step
Nsim = ceil(65 / dt); % Total simulation time + Planning time
mode = 1; %% mode 1= zigzag, mode 2= circling
param.V_c = 0; % Ocean current velocity

% mode = test_param1(iii);
% param.V_c = test_param2(jjj);
%%
tic;
param.Nsim = Nsim;
param.Nx = 12; % dimension x
param.Nu = 2; % dimension of u
param.Nt = 100; % horizon length
param.dt = dt; % step time
param.xmax = inf(param.Nx, 1); % state constraints
param.xmin = -param.xmax;
param.mpc_cal_t = 100; % computation limit of NMPC

I = [17.1572 0 4.4 0 11 0; % Mass Matrix
    0 37.07 0 -11 0 -11;
    4.4 0 45.265 0 11 0;
    0 -11 0 60.2815 0 0;
    11 0 11 0 137.5 0;
    0 -11 0 0 0 110];

param.umax = ones(param.Nu, 1) * 115; % input constraints
param.umin = -ones(param.Nu, 1) * 65;
param.I = I; % body inertial tensor
%% generate reference trajectory
q0_ref = [1;0;0;0];
p0_ref = [0;0;0];
w0_ref = [0;0;0.1] ; % desired yaw velocity 0.1 rad/s
v0_ref = [0.5;0;0] ; % desired surge velocity 0.5 m/s
x0_ref = [q0_ref;p0_ref];
xid_ref = [w0_ref;v0_ref];

X = eye(4);

X_ref = x0_ref';
xi_ref = xid_ref';
for i = 1:Nsim
    xid_ref_rt = xid_ref';

    if mode == 1 % zigzag manuever
        xid_ref_rt(3) = 0.1*cos(i/100); 
    end
    Xi = [skew(xid_ref_rt(1:3)), xid_ref_rt(4:6)';...
        [0,0,0,0]];
    X = X * expm(Xi * dt);
    X_ref = [X_ref; [rotm2quat(X(1:3,1:3)), X(1:3,4)']];
    xi_ref = [xi_ref; xid_ref_rt];
end
param.X_ref = X_ref;
param.xi_ref = xi_ref;
%%
Logger1 = struct('X',zeros(13, Nsim - 11),'Err',zeros(6, Nsim - 12),'U',zeros(6, Nsim - 12),'Err2', zeros(Nsim - 12, 4), 'taken_time', zeros(1));
Logger2 = Logger1;
Logger3 = Logger1;

parfor k = 1:Ns
    % Setting initial states
    q0 = eul2quat([0, 0, (rand-0.5)*0.5*pi], 'XYZ')';
    p0 = zeros(3,1);
    p0(1:2) = p0(1:2) + (rand(2,1)-1)*5;
    p0(2) = p0(2) + 2.5;
    w0 = (rand(3,1) - 0.5) * 0;
    v0 = (rand(3,1) - 0.5) * 0;
    v0(1) = 0.5;

    Logger1(k) = sim_mpc(q0, p0, w0, v0, dt, Nsim, 1, param); % Proposed MPC
    Logger2(k) = sim_mpc(q0, p0, w0, v0, dt, Nsim, 2, param); % NMPC
    Logger3(k) = sim_mpc(q0, p0, w0, v0, dt, Nsim, 3, param); % NMPC-simple
end
plot_bool = 0;
if plot_bool == 1
%% Error plots for all trajectories

axis_font_size = 12;
lw_1 = 0.1;

x_lim = [0, 60];
figure(1)
subplot(3,2,1)
for k = 1:Ns
    plot(0:dt:(length(Logger1(k).Err2)-1) * dt, Logger1(k).Err2(:,1),'LineWidth',lw_1)
    hold on
end
title("Proposed MPC - Orientation error", "Interpreter","latex")
box on
grid on
xlim(x_lim)
set(gca,'FontSize',axis_font_size)

subplot(3,2,3)
for k = 1:Ns
    plot(0:dt:(length(Logger2(k).Err2)-1) * dt, Logger2(k).Err2(:,1),'LineWidth',lw_1)
    hold on
end
title("NMPC - Orientation error", "Interpreter","latex")
ylabel("$\|Log(R_d^{-1}R)\|$", "Interpreter","latex")
box on
grid on
xlim(x_lim)
set(gca,'FontSize',axis_font_size)

subplot(3,2,5)
for k = 1:Ns
    plot(0:dt:(length(Logger3(k).Err2)-1) * dt, Logger3(k).Err2(:,1),'LineWidth',lw_1)
    hold on
end
title("NMPC-simple - Orientation error", "Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
box on
grid on
xlim(x_lim)
set(gca,'FontSize',axis_font_size)

subplot(3,2,2)
for k = 1:Ns
    plot(0:dt:(length(Logger1(k).Err2)-1) * dt, vecnorm(Logger1(k).Err2(:,2:4),2,2),'LineWidth',lw_1)

    hold on
end
title("Proposed MPC - Position error", "Interpreter","latex")
box on
grid on
xlim(x_lim)
set(gca,'FontSize',axis_font_size)

subplot(3,2,4)
for k = 1:Ns
    plot(0:dt:(length(Logger2(k).Err2)-1) * dt, vecnorm(Logger2(k).Err2(:,2:4),2,2),'LineWidth',lw_1)
    hold on
end
title("NMPC - Position error", "Interpreter","latex")
ylabel("$\|p-p_d\|$", "Interpreter","latex")
box on
grid on
xlim(x_lim)
set(gca,'FontSize',axis_font_size)

subplot(3,2,6)
for k = 1:Ns
    plot(0:dt:(length(Logger2(k).Err2)-1) * dt, vecnorm(Logger3(k).Err2(:,2:4),2,2),'LineWidth',lw_1)
    hold on
end
title("NMPC-simple - Position error", "Interpreter","latex")
xlabel("Time (s)", "Interpreter","latex")
box on
grid on
xlim(x_lim)
set(gca,'FontSize',axis_font_size)

%% Trajectory plot

figure(2)
loggers = {Logger1, Logger2, Logger3};
k = 9; % i-th trajectory
for id = 1:3
    lw = 2;
    logger = loggers{id};
    if id == 1
        color = "r";
    elseif id == 2
        color = "g";
    elseif id == 3
        color = "b";
    end
       
    plot(logger(k).X(5, 1:1200), logger(k).X(6, 1:1200), color+"-.", "LineWidth", lw)
    hold on
end
plot(X_ref(1:1200,5), X_ref(1:1200,6), '-','color','k' ,"LineWidth",3)
daspect([1,1,1])
box on
grid on
legend({"Proposed MPC", "NMPC", "NMPC-simple", "Rererence"}, "interpreter", "latex")
xlabel("x", "interpreter", "latex")
ylabel('y', "interpreter", "latex")
zlabel('z', "interpreter", "latex")

%%
end
