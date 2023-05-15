function [x, q] = rk4_(h, f, x, Q, q)
% [xplus, quad] = rk4_int_(h, f, x, [Q], [q])
%
% Low-level interface for one RK4 step. f must be a function of only x.
narginchk(3, 5)
if nargin() == 4
    error('Must provide Q and q for quadrature!');
end
k1 = f(x);
k2 = f(x + k1*h/2);
k3 = f(x + k2*h/2);
k4 = f(x + k3*h);
if nargin() >= 4
    k1_q = Q(x);
    k2_q = Q(x + k1*h/2);
    k3_q = Q(x + k2*h/2);
    k4_q = Q(x + k3*h);
    q = q + (k1_q + 2*k2_q + 2*k3_q + k4_q)*h/6;
else
    q = 0;
end
x = x + (k1 + 2*k2 + 2*k3 + k4)*h/6;
end%function

