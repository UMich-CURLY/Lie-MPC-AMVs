function [A, bmin, bmax] = eMPCConstraints(xid, p0, x0, param)
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
D = diag([-54.3823, -246.0496, -45.265 *(1+10*abs(x0(10))),-77.5544, 0, -546.4805]);
I_inv = I^(-1);

for k = 1:Nt
    D_t = D;
    D_t(3) = D_t(3) - 45.265 *(10*abs(xid(k,3)));

    Irb =  [ 55.0000         0         0         0  -11.0000         0;
         0   55.0000         0   11.0000         0   11.0000;
         0         0   55.0000         0  -11.0000         0;
         0   11.0000         0   14.6643         0    4.4000;
  -11.0000         0  -11.0000         0   22.5500         0;
         0   11.0000         0    4.4000         0   18.1500];
    M11 = Irb(1:3, 1:3)  ;
    M12 = Irb(1:3, 4:6) ;
    M21 = Irb(4:6, 1:3)  ;
    M22 = Irb(4:6, 4:6)  ;

    xi_w = xid(k,1:3)';
    xi_v = xid(k,4:6)';
    G = zeros(6,6);
    G(1:3,1:3) = skew(xi_w)*M22 + skew(xi_v)*M12;
    G(1:3,4:6) = skew(xi_v)*M11 + skew(xi_w)*M21;
    G(4:6,1:3) = skew(xi_w)*M12;
    G(4:6,4:6) = skew(xi_w)*M11;

    Crb = zeros(6,6);
    Crb(1:3,1:3) = skew(M21*xi_v) + skew(M22*xi_w);
    Crb(1:3,4:6) = skew(M11*xi_v) + skew(M12*xi_w);
    Crb(4:6,1:3) = skew(M11*xi_v) + skew(M12*xi_w);
    
    Iam = diag([5.2815, 82.5, 55, 2.4929, 14.52, 27.115]);
    Mam11 = Iam(1:3,1:3)  ;
    Mam22 = Iam(4:6,4:6)  ;

    Cam = zeros(6,6);
    Cam(1:3,1:3) = skew(Mam22*xi_w);
    Cam(1:3,4:6) = skew(Mam11*xi_v);
    Cam(3,4) = 0;
    Cam(3,5) = 0;
    Cam(4:6,1:3) = skew(Mam11*xi_v);
    Cam(4,3) = 0;
    Cam(5,3) = 0;

    G2 = zeros(6,6);  
    G2(1:3,1:3) = skew(xi_w)*Mam22;
    G2(1:3,4:6) = skew(xi_v)*Mam11 ;
    G2(4:6,4:6) = skew(xi_w)*Mam11;
    G2(3,4) = 0;
    G2(3,5) = 0;
    G2(4,6) = 0;
    G2(5,6) = 0;

    Cam= -Cam;
    Crb = -Crb;

    H = I_inv *(Crb + Cam + G + G2 + D_t);
    
    Ac = [-adjoint_(xid(k, :)),   -eye(6);...
        zeros(6), H];
    Bc = [zeros(6,2);...
        I_inv * [0, 0; 0, 0; 0.395, -0.395; 1, 1; 0, 0; 0, 0]];

    b =  (G + G2) * x0(8:13); 
    b(3) = b(3) + 45.265 *(10*abs(xid(k,3)))*xid(k,3);
    
    b = -I_inv * b;
    hc = [xid(k,:)'; b];

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