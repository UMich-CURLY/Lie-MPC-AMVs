function [xhat, status] = nlfilter(varargin)
% `[xhat, status] = nlfilter(h, xhatm, y, Rinv, Pinv, [xlb], [xub])`
%
% Solves the nonlinear filtering problem
%
% > $\min_x (x - \hat{x}^-)'P^{-1}(x - \hat{x}^-) + (y - h(x))'R^{-1}(y - h(x))$
%
% This is essentially a zero-step MHE problem with quadratic prior and cost.
%
% Arguments are as follows:
%
% * `h` : Casadi Function
%
%   > Measurement function giving $y = h(x)$.
%
% * `xhatm` : Column Vector
%
%   > Numerical estimate for $\hat{x}^-$, i.e., $\hat{x}(k | k - 1)$.
%
% * `y` : Column Vector
%
%   > Numerical value for the current measurement $y(k)$.
%
% * `Rinv` : Matrix
% * `Pinv` : Matrix
%
%   > Penalty matrices $R^{-1}$ and $P^{-1}$. Note that they should already be
%   inverted when they are passed as arguments.
%
% * `xlb` : Column Vector
% * `xub` : Column Vector
%
%   > Bounds to enforce on $x$ in the optimization problem. Useful if $h$ is
%   undefined for certain values of $x$ and you need to restrict the domain.
%
% The two outputs are the column vector `xhat`, which is the optimal estimate
% $\hat{x}(k | k)$ and a string `status` that gives the status of the
% optimization problem.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('h', 'required', 'casadifunc');
    parser.add('xhatm', 'required', {'numeric', 'col'});
    parser.add('y', 'required', {'numeric', 'col'});
    parser.add('Rinv', 'required', 'numeric');
    parser.add('Pinv', 'required', 'numeric');
    parser.add('xlb', -inf(), {'numeric', 'col'});
    parser.add('xub', inf(), {'numeric', 'col'});
end
args = parser.parse(varargin{:});

% Make variable and objective function.
x = casadi.SX.sym('x', length(args.xhatm));
xbar = args.xhatm;
dx = x - xbar;
dy = args.y - args.h(x);
f = dx'*args.Pinv*dx + dy'*args.Rinv*dy;

% Choose bounds.
xlb = args.xlb.*ones(size(xbar));
xub = args.xub.*ones(size(xbar));

% Build nlp and solver.
nlp = struct('x', x, 'f', f);
solver = casadi.nlpsol('nlfilter', 'ipopt', nlp, ...
                       struct('print_time', false(), ...
                              'ipopt', struct('print_level', 0)));
sol = solver('x0', xbar, 'lbx', xlb, 'ubx', xub);
xhat = full(sol.x);

stats = solver.stats();
status = stats.return_status;

end%function

