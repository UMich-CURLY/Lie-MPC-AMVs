function [A, bmin, bmax] = eMPC_nodamp_Constraints(xid, p0, x0, param)
% min x'Qx + u'Ru
% xkp1 = Ak xk + Bk uk + bk, k = 0,1, ..., N-1
% x0 = x_init,
% umin < uk < umax, k = 0, 1, ..., N-1

I = param.I;
dt = param.dt;

Nx = param.Nx;
Nu = param.Nu;
Nt = param.Nt;
A = zeros(Nx * (Nt + 1) + Nu * Nt, Nt * (Nx + Nu) + Nx);
Noff =  Nx * (Nt + 1);
bmin =zeros(size(A,1), 1);

for k = 1:Nt
    % xi_bar = I * xid(k, :)';
    xi_bar = I * x0(8:13);
    G = zeros(6,6);
    G(1:3,1:3) = skew(xi_bar(1:3));
    G(1:3,4:6) = skew(xi_bar(4:6));
    G(4:6,1:3) = skew(xi_bar(4:6));

    H = I^-1 *(coadjoint(x0(8:13)) * I + G);

    Ac = [-adjoint_(xid(k, :)),   -eye(6);...
        zeros(6), H];
    Bc = [zeros(6,2);...
        I^-1 * [0, 0; 0, 0; 0.395, -0.395; 1, 1; 0, 0; 0, 0]];

    b = -I^(-1) * G * x0(8:13);
    hc = [xid(k,:)'; b];

    %     Ad = expm(Ac * dt);
    %     Bd = (eye(12) * dt + Ac * dt^2 / 2 + Ac^2 * dt^3 / 6) * Bc;
    %     hd = (eye(12) * dt + Ac * dt^2 / 2 + Ac^2 * dt^3 / 6) * hc;
    Ad = eye(param.Nx) + Ac * dt;
    Bd = dt * Bc;
    hd = dt * hc;

    A((k-1)*Nx+1:k*Nx, (k-1)*Nx+1:k*Nx) = -Ad;
    A((k-1)*Nx+1:k*Nx, k*Nx+1:(k+1)*Nx) = eye(Nx);
    A((k-1)*Nx+1:k*Nx, (Noff+(k-1)*Nu+1):(Noff+k*Nu)) = -Bd;
    bmin(Nx * (k-1)+1:Nx * k) = hd;
end
A(Nt*Nx+1:(Nt+1)*Nx, 1:Nx) = eye(Nx);
for k = 1:Nt
    A(Noff+1+(k-1)*Nu:Noff+k*Nu,  Noff+1+(k-1)*Nu:Noff+k*Nu) = eye(Nu);
end
bmin(Nt*Nx+1:(Nt+1)*Nx) = p0(:);
bmax = bmin;
bmax((Nt+1)*Nx+1:end) = repmat(param.umax(:), [Nt, 1]);
bmin((Nt+1)*Nx+1:end) = repmat(param.umin(:), [Nt, 1]);
end