function [Pmp1, xhatmp1, P, xhat] = ekf(varargin)
% `[Pm, xhatm, P, xhat] = ekf(f, h, x, u, w, y, Pm, Q, R, [projectionfunc])`
%
% Updates the prior distribution $P^-$ using the Extended Kalman filter.
%
% `f` and `h` should be Casadi functions. `f` must be discrete-time. P, Q, and R
% are the prior, state disturbance, and measurement noise covariances (in
% particular, the input `Pm` refers to $P(k | k - 1)$).
%
% Note that `f` must be $f(x,u,w)$ and `h` must be $h(x)$.
%
% The value of `x` that should be fed is $\hat{x}(k | k-1)$, and the value of
% `P` should be P(k | k-1). xhat will be updated to xhat(k | k) and then
% advanced to $\hat{x}(k+1 | k)$, while `P` will be updated to $P(k | k)$
% and then advanced to $P(k+1 | k)$. The return values are a list as follows
%
% > $[P(k+1 | k), \hat{x}(k+1 | k), P(k | k), \hat{x}(k | k)]$
%    
% Depending on your specific application, you will only be interested in
% some of these values.
%
% The optimal argument `projectionfunc` is a function handle that is called to
% project $x$ back into the feasible space after the correction step. It should
% take $x$ as its only argument and return the projected version of $x$. For
% example, if $x$ is nonnegative, using
%
% > `projectionfunc = @(x) max(x, 0)`
%
% would project x back into the nonnegative orthant before advancing.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('f', 'required', 'casadifunc');
    parser.add('h', 'required', 'casadifunc');
    parser.add('x', 'required', {'numeric', 'col'});
    parser.add('u', 'required', {'numeric', 'col'});
    parser.add('w', 'required', {'numeric', 'col'});
    parser.add('y', 'required', {'numeric', 'col'});
    parser.add('Pm', 'required', 'numeric');
    parser.add('Q', 'required', 'numeric');
    parser.add('R', 'required', 'numeric');
    parser.add('projectionfunc', []);
end
args = parser.parse(varargin{:});

% Get jacobians.
f = args.f;
f_jacx = mpctools.jacobian(f);
f_jacw = mpctools.jacobian(f, 3);

h = args.h;
h_jacx = mpctools.jacobian(h);
    
% Get linearization of measurement.
x = args.x;
C = full(h_jacx(x));
yhat = full(h(x));

% Advance from x(k | k-1) to x(k | k).
xhatm = args.x; % This is xhat(k | k-1)    
Pm = args.Pm; % This is P(k | k-1)

L = ((C*Pm*C' + args.R)\(C*Pm))';

xhat = xhatm + L*(args.y - yhat); % This is xhat(k | k) 
P = (eye(size(Pm)) - L*C)*Pm; % This is P(k | k)

% Perform correction.
if ~isempty(args.projectionfunc)
    xhat = args.projectionfunc(xhat);
end

% Now linearize the model at xhat.
w = zeros(size(args.w));
u = args.u;
A = full(f_jacx(xhat, u, w));
G = full(f_jacw(xhat, u, w));

% Advance.
Pmp1 = A*P*A' + G*args.Q*G'; % This is P(k+1 | k)
xhatmp1 = full(f(xhat, u, w)); % This is xhat(k+1 | k)    

end%function

