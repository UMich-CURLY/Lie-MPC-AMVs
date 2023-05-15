function [A, bmin, bmax] = MPCConstraints(Ad, Bd, b, x0, param)
% min x'Qx + u'Ru
% xkp1 = Ak xk + Bk uk + bk, k = 0,1, ..., N-1
% x0 = x_init,
% umin < uk < umax, k = 0, 1, ..., N-1
Nx = param.Nx;
Nu = param.Nu;
Nt = param.Nt;
A = zeros(Nx * (Nt + 1) + Nu * Nt, Nt * (Nx + Nu) + Nx);
Noff =  Nx * (Nt + 1);
for k = 1:Nt
   A((k-1)*Nx+1:k*Nx, (k-1)*Nx+1:k*Nx) = -Ad;
   A((k-1)*Nx+1:k*Nx, k*Nx+1:(k+1)*Nx) = eye(Nx);
   A((k-1)*Nx+1:k*Nx, (Noff+(k-1)*Nu+1):(Noff+k*Nu)) = -Bd;
end
A(Nt*Nx+1:(Nt+1)*Nx, 1:Nx) = eye(Nx);
for k = 1:Nt
    A(Noff+1+(k-1)*Nu:Noff+k*Nu,  Noff+1+(k-1)*Nu:Noff+k*Nu) = eye(Nu);
end
bmin =zeros(size(A,1), 1);
bmin(1:Nt * Nx) = repmat(b, [Nt, 1]);
bmin(Nt*Nx+1:(Nt+1)*Nx) = x0(:);
bmax = bmin;
bmax((Nt+1)*Nx+1:end) = repmat(param.umax(:), [Nt, 1]);
bmin((Nt+1)*Nx+1:end) = repmat(param.umin(:), [Nt, 1]);

% velocity constraints. 
Av = zeros( (Nx/2) * Nt, Nt * (Nx + Nu) + Nx);
for k = 2:Nt+1
    Av((k-2)*6+1:6*(k-1), Nx * (k - 1) + 6+1:Nx * (k)) = eye(6);
end
bvmin = -ones(6,1) * 1.5;
bvmin = repmat(bvmin, [Nt, 1]);
bvmax = -bvmin;

% A = [A; Av];
% bmin = [bmin; bvmin];
% bmax = [bmax; bvmax];
end