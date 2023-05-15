function [X,Xuu,Xu] = forceSurgeDamping(flag,u_r,m,S,L,T1,rho,u_max,thrust_max)
% [X,Xuu,Xu] = forceSurgeDamping(flag,u_r,m,S,L,T1,rho,u_max,thrust_max) 
% computes the hydrodynamic damper in surge for low- and high-speed marine 
% craft given by
%
%  (m - Xudot) d/dt nu_r + d(nu_r) * nu_r = thrust 
%
% where nonlinear damping d(nu_r) * nu_r = -X. The linear and quadratic
% damping terms are combined using a blending function 
% sigma = exp( -alpha * u_r^2 ) with alpha = 5, such that
%
%    X = sigma * Xu * u_r    +  (1 - sigma) * Xuu * abs(u_r) * u_r
%           (linear damping)             (quadratic damping)  
%
% The blending function ensures that "double counting" is avoided for
% ships/vehicles that operate at low speed (dominated by linear damping) and 
% high speed (dominated by quadratic damping). This is verfified by
%
% low speed  (sigma = 1): linear damping dominates, X = Xu * u_r
% high speed (sigma = 0): quadratic damping dominates, X = Xuu * abs(u_r)*u_r 
%
% Outputs:
%  X:   total damping force in surge (N)
%  Xu:  OPTIONAL linear damping coefficient
%  Xuu: OPTIONAL quadratic damping coefficient
%
% Inputs:
%  flag:  1 - plot/inspect the linear and quadratic damping terms
%         0 - supress plotting, numerical ouputs only
%  u_r:   relative surge velocity where u_r = u - u_c (m/s) 
%  m:     mass of the ship/vehicle (kg)
%  S:     wetted surface, area of the outer surface in contact with the water (m^2)
%  L:     length of ship/vehicle (m)
%  T1:    time constant in surge used to estimate the linear damping coeff. (s)
%  rho:   density of water (kg/m3)
%  u_max: maximum surge velocity of the ship/vehicle, no currents/waves (m/s)
%  thrust_max: OPTIONAL parameter specyfying the maximum thrust in surge 
%         corresponding to u_max. If specified, the quadratic damping
%         coefficient Xuu is calibrated. If the argument is omitted, Xuu
%         is computed using the ITTC resistance curve 
%
% Examples:
%  Quadratic damping based on the ITTC resistance curve:
%   X = forceSurgeDamping(flag,u_r,m,S,L,T1,rho,u_max)
%  Quadratic damping obtained from the steady-state condition 
%  -Xuu * u_max^2 = thrust_max using the optional input thrust_max:
%   X = forceSurgeDamping(flag,u_r,m,S,L,T1,rho,u_max,thrust_max)
%
% Author:    Thor I. Fossen
% Date:      17 Dec 2021 

% linear damping coefficient
Xudot = -addedMassSurge(m,L,rho);
Xu = -(m - Xudot) / T1;      

% quadratic damping coefficient
if (nargin == 9)  
    
    % steady state condition: -Xuu * u_max^2 = thrust_max
    Xuu = - thrust_max / u_max^2;
    
else

    nu_kin = 1e-6;                                % kinematic visosity
    k = 0.1;                                      % correction factor
    eps = 1e-10;                                
    Rn = (L / nu_kin) * abs(u_r);                 % Reynolds number
    Cf = 0.075 / (log10(Rn + eps) -  2)^2;        % ITTC resistance curve
    
    Xuu = -0.5 * rho * S * (1-k) * Cf;
    
end

% blending function 
alpha = 5;
sigma = exp( -alpha * u_r^2 );

% damping force in surge
X = sigma * Xu * u_r + (1 - sigma) * Xuu * abs(u_r) * u_r;

% Plot damping terms when flag = 1
if (flag == 1)
    u = 0:0.1:u_max+1;  
    sigma = exp( -alpha * u.^2 );
    T_max = -Xuu * u_max^2;
    d_lin = -sigma .* Xu; 
    d_quad = - (1 - sigma) .* Xuu .* abs(u);
    figure(gcf)
    subplot(211)
    plot(u,(d_lin+d_quad).*u,'-',u_max,T_max,'*','linewidth',2), grid
    legend('Linear + quadratic damping force (N)','Max speed/thrust');
    title('X = sigma * Xu * u + (1 - sigma) * Xuu * abs(u) * u')
    subplot(212)    
    u = 0:0.1:3;  
    sigma = exp( -alpha * u.^2 );
    plot(u,-sigma * Xu,u,-(1 - sigma) .* Xuu .* abs(u),'linewidth',2), grid;
    legend('Linear damping coefficient: -Xu','Quadratic damping term: -Xuu * |u|');
end

end