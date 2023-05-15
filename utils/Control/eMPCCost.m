function [M, q, solver] = eMPCCost(Q, R, P, xid, param, solver)
% min x'Qx + (xi-xid)'Q(xi-xid) + u'Ru
% xkp1 = Ak xk + Bk uk + bk, k = 0,1, ..., N-1
% x0 = x_init,
% umin < uk < umax, k = 0, 1, ..., N-1

q = zeros((param.Nt+1)*param.Nx + param.Nt * param.Nu,1);
M = zeros((param.Nt+1)*param.Nx + param.Nt * param.Nu,(param.Nt+1)*param.Nx + param.Nt * param.Nu );
Nx = param.Nx;

for k = 1:param.Nt-1
    %% consider the transport map
    C = eye(12);
    C(7:12, 1:6) = -adjoint_(xid(k,:));
    M((k-1)*Nx+1:(k)*Nx, (k-1)*Nx+1:(k)*Nx) = C'*Q'*C;
    b = zeros(12,1);
    b(7:12) = xid(k,:)';
    q((k-1)*Nx+1:(k)*Nx, 1) = -C' * Q * b;
end
k = param.Nt;
C = eye(12);
C(7:12, 1:6) = -adjoint_(xid(k,:));
b = zeros(12,1);
b(7:12) = xid(k,:)';
b = C' * P * b;
M((k-1)*Nx+1:(k)*Nx,(k-1)*Nx+1:(k)*Nx) = C'*P*C;
q((k-1)*Nx+1:(k)*Nx) = -b; % param.Ad *
Nu =param.Nu;
for j = 1:param.Nt
    M((k+1)*Nx+(j-1)*Nu+1:(k+1)*Nx+(j)*Nu, (k+1)*Nx+(j-1)*Nu+1:(k+1)*Nx+(j)*Nu) = R;
end
end