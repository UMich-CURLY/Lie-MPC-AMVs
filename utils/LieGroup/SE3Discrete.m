function [Ad, Bd, b] = SE3Discrete(x,param)
% [φ, θ, ψ, x, y, z, ω1, ω2, ω3, vx, vy, vz]
r = [x(1);x(2);x(3)];
p = [x(4);x(5);x(6)];
w = [x(7);x(8);x(9)];
v = [x(10);x(11);x(12)];


Ad = [eye(6), eye(6) * param.dt;...
     zeros(6), eye(6)];
% ZOH? Replace later. 
 Bd = [J^(-1) * dt^2;
     J^(-1) * dt];
 b = [1/2 * dt^2]

end