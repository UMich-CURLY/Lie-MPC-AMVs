function dxdt = SE3Dyn(t,x,u,I)
% [φ, θ, ψ, x, y, z, ω1, ω2, ω3, vx, vy, vz]
q = [x(1);x(2);x(3);x(4)];
p = [x(5);x(6);x(7)];
w = [x(8);x(9);x(10)];
v = [x(11);x(12);x(13)];

Omega = [0, -w(1), -w(2), -w(3);...
         w(1), 0, w(3), -w(2);...
         w(2), -w(3), 0, w(1);...
         w(3), w(2), -w(1), 0];
dQuat = 1 / 2 * Omega * q;
R = quat2rotm(q');
dx = [dQuat; R * v];
d2x = I^(-1) *  (coadjoint([w;v]) * I * [w;v] + u + 0 * I *  [0;0;0;R'*[0;0;-9.81]]);
dxdt = [dx; d2x];
end