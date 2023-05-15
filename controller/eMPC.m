function [u, solver] = eMPC(x0, xd, xid, param, solver)
Xd = [quat2rotm(xd(1:4)), xd(5:7)';...
    0,0,0,1];
X00 = [quat2rotm(x0(1:4)'), x0(5:7);...
    0,0,0,1];

% The inverse of Left error defined in paper is equivalent. 
X00 =  X00^(-1) * Xd; 
X0 = logm(X00);

p0 = [X0(3,2); X0(1,3); X0(2,1); X0(1:3,4);x0(8:end)]; %% error

[A, bmin, bmax] = eMPCConstraints(xid, p0, x0, param);
[M, q, solver] = eMPCCost(param.Q, param.R, param.P, xid, param, solver);
    
solver.solver = osqp;
solver.solver.setup(M, q, A, bmin, bmax,'verbose', 0);

sol = solver.solver.solve;
u = sol.x((param.Nt+1)*param.Nx+1:(param.Nt+1)*param.Nx+param.Nu);
end
