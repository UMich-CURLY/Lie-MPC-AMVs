function [M, q] = gMPCCost(Q, R, P, xid, param)
% min x'Qx + (xi-xid)'Q(xi-xid) + u'Ru
% xkp1 = Ak xk + Bk uk + bk, k = 0,1, ..., N-1
% x0 = x_init,
% umin < uk < umax, k = 0, 1, ..., N-1 
M = Q * 0;
q = zeros((param.Nt+1)*param.Nx + param.Nt * param.Nu,1);
for k = 1:param.Nt-1
    M = blkdiag(M, Q);
    w = diag(Q);
    % q((k-1)*param.Nx+7:(k)*param.Nx) = - w(7:end) .* (xid(k,:)'); % param.Ad *
end
M = blkdiag(M, P);
w = diag(P);
k = param.Nt;
% q((k-1)*param.Nx+7:(k)*param.Nx) = - w(7:end) .*  (xid(k,:)'); % param.Ad * 
offset = param.Nx * (param.Nt + 1);
for k = 1:param.Nt
    M = blkdiag(M, R);
    ud = -coadjoint(xid(k,:)) * param.I * xid(k,:)';
    ud = ud([4,5,6,1,2,3]);
    q(offset + (k-1)*param.Nu+1:offset+(k)*param.Nu) = - R * ud; % param.Ad *
end
end