function [first, second] = rk4(varargin)
% `x = rk4(f, x0, [par={}], [Delta=1], [M=1])`
%
% `[q, x] = rk4(..., 'quad', Q)`
%
% Does `M` steps of explicit rk4 integration with a total time of `Delta`.
%
% * `f` : Function
%
%   > ODE right-hand side function.
%
% * `x0` : Vector
%
%   > Initial condition.
%
% * `par` : Cell Array
%
%   > Parameters (extra arguments) to `f` (called as `f(x0, par{:})`).
%
% * `Delta` : Positive Scalar
%
%   > Timestep to use for integration.
%
% * `M` : Positive Integer
%
%   > Number of steps to take. Note that the total time is always `Delta`,
%   and so the duration of each step is $\Delta/M$.
%
% * `quad` : Function
%
%   > Quadrature function. If provided, must take the same arguments as f.
%   Note that this makes the first return value equal to the quadrature rather
%   than the state at the next time point.
%
% Returned values are `x`, the approximate value of `x` after `Delta` time
% units; and `q`, the approximate integral of `Q(x, ...)` over the timestep.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('f');
    parser.add('x0');
    parser.add('par', {}, {'cell', 'row'});
    parser.add('Delta', 1, {'scalar', 'pos'});
    parser.add('M', 1, {'scalar', 'pos', 'int'});
    parser.add('quad', []);
end
args = parser.parse(varargin{:});
needquad = ~isempty(args.quad);

% Perform RK4 steps.
h = args.Delta/args.M;
par = args.par;
f = @(x) args.f(x, par{:});
if needquad
    Q = @(x) args.Q(x, par{:});
end
x = args.x0;
q = 0;
for j = 1:args.M
    if needquad
        [x, q] = mpctools.rk4_(h, f, x, Q, q);
    else
        x = mpctools.rk4_(h, f, x);
    end
end

% Decide output arguments.
if needquad
    first = q;
    second = x;
else
    first = x;
    second = [];
end

end%function
