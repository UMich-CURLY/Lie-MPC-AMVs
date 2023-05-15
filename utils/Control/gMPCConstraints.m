function [A, bmin, bmax] = gMPCConstraints(xid, p0, R0, param)
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
        dRR = expm(dt * skew(xid(1:3)));
        Ad = [eye(3), dt * eye(3),  zeros(3), zeros(3);...
            zeros(3),      eye(3),  zeros(3), zeros(3);...
            zeros(3),     zeros(3),     dRR', dt * dRR';...
            zeros(3),     zeros(3), zeros(3), I(1:3,1:3)^(-1) * dRR' * (I(1:3,1:3) + dt * skew(I(1:3,1:3) * xid(1:3)'))];
        Bd = [zeros(3), zeros(3);...
            dt * I(4:6,4:6)^(-1) * R0, zeros(3);...
            zeros(3), zeros(3);...
            zeros(3), dt *  I(1:3,1:3)^(-1)];
        A((k-1)*Nx+1:k*Nx, (k-1)*Nx+1:k*Nx) = -Ad;
        A((k-1)*Nx+1:k*Nx, k*Nx+1:(k+1)*Nx) = eye(Nx);
        A((k-1)*Nx+1:k*Nx, (Noff+(k-1)*Nu+1):(Noff+k*Nu)) = -Bd;
        bd = zeros(Nx, 1);
        bmin(Nx * (k-1)+1:Nx * k) = bd;
    end
    A(Nt*Nx+1:(Nt+1)*Nx, 1:Nx) = eye(Nx);
    for k = 1:Nt
        A(Noff+1+(k-1)*Nu:Noff+k*Nu,  Noff+1+(k-1)*Nu:Noff+k*Nu) = eye(Nu);
    end
    % bmin(1:Nt * Nx) = repmat(b, [Nt, 1]);
    bmin(Nt*Nx+1:(Nt+1)*Nx) = p0(:);
    bmax = bmin;
    bmax((Nt+1)*Nx+1:end) = repmat(param.umax(:), [Nt, 1]);
    bmin((Nt+1)*Nx+1:end) = repmat(param.umin(:), [Nt, 1]);
    
    %     % velocity constraints.
    %     Av = zeros( (Nx/2) * Nt, Nt * (Nx + Nu) + Nx);
    %     for k = 2:Nt+1
    %         Av((k-2)*6+1:6*(k-1), Nx * (k - 1) + 6+1:Nx * (k)) = eye(6);
    %     end
    %     bvmin = -ones(6,1) * 1.5;
    %     bvmin = repmat(bvmin, [Nt, 1]);
    %     bvmax = -bvmin;
    %
    % A = [A; Av];
    % bmin = [bmin; bvmin];
    % bmax = [bmax; bvmax];
end