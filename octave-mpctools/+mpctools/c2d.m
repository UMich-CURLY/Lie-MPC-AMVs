function varargout = c2d(varargin)
% `[...] = c2d(Delta, a, b, [q, r], [m], [g, h], [quad], [Nquad], [return])`
%
% Discretization with continuous objective.
%
% Converts from continuous-time objective
%
% > $l(x,u) = \int_0^\Delta x'qx + 2x'mu + u'ru + g'x + h'u \; dt \qquad$
% > $dx/dt = ax + bu$
%    
% to the equivalent (assuming `u` is constant on $[0, \Delta]$).
%
% > $L(x,u) = x'Qx + 2x'Mu + u'Qu + G'x + H'u \qquad$
% > $x^+ = Ax + Bu$
%
% in discrete time.
%
% Note that `q` can be given if and only if `r` is given (similarly `g` and
% `h`).
%
% Optional argument `quad` decides whether to use approximate quadrature to
% compute the Q, R, and M matrices (default `false`). `Nquad` is the number of
% steps to use in the quadrature (default 100). Quadrature may be necessary
% when `Delta*a` has eigenvalues with real part less than $-25$ or greater than
% $25$.
%
% Argument `return` decides what is returned. The default value is 'struct',
% which returns a struct with fields "A" and "B", as well as fields "Q", "R",
% and "M" if arguments `q` and `r` were given. `return` can also be a string
% consisting of "ABQRMGH", at which point, the individual matrices given in the
% string will be returned.
%
% Reference: C. Van Loan, 1978, "Computing integrals involving the matrix
% exponential".
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('Delta', 'required', {'scalar', 'pos'});
    parser.add('a', 'required', 'numeric');
    parser.add('b', 'required', 'numeric');
    parser.add('q', [], 'numeric');
    parser.add('r', [], 'numeric');
    parser.add('m', [], 'numeric');
    parser.add('g', [], 'numeric');
    parser.add('h', [], 'numeric');
    parser.add('quad', false(), 'bool');
    parser.add('Nquad', 100, {'scalar', 'int'});
    parser.add('return', 'struct', 'str');
end
args = parser.parse(varargin{:});

% Check sizes.
Nx = size(args.a, 1);
if size(args.a, 2) ~= Nx
    error('a must be square!');
end
Nu = size(args.b, 2);
if size(args.b, 1) ~= Nx
    error('b must be Nx by Nu!');
end

% Check if objective.
needobjective = ~isempty(args.q) || ~isempty(args.r) || ~isempty(args.m);
if needobjective
    if isempty(args.q) || isempty(args.r)
        error('Both q and r must be provided if either is provided!');
    end
    if isempty(args.m)
        args.m = zeros(Nx, Nu);
    end
end

% Use matrix exponential formulas with or without objective terms.
mats = struct();
if needobjective && ~args.quad
    [i, E] = getindices([Nx, Nu, Nx, Nu]);
    E(i{1},i{1}) = -args.a';
    E(i{1},i{3}) = args.q;
    E(i{1},i{4}) = args.m;
    E(i{2},i{1}) = -args.b';
    E(i{2},i{3}) = args.m';
    E(i{2},i{4}) = args.r;
    E(i{3},i{3}) = args.a;
    E(i{3},i{4}) = args.b;
    
    % Exponentiate and grab pieces.
    [E, nsquare] = expmnosquare(E*args.Delta);
    A = E(i{3},i{3});
    B = E(i{3},i{4});
    Q = A'*E(i{1},i{3});
    R = B'*E(i{1},i{4}) + E(i{2},i{4});
    M = A'*E(i{1},i{4});
    
    % Use recursive formulas for squaring.
    for k = 1:nsquare
        R = 2*R + B'*M + M'*B + B'*Q*B;
        M = M + A'*(Q*B + M);
        Q = Q + A'*Q*A;
        B = B + A*B;
        A = A^2;
    end
    
    % Package up matrices, adding extra part for R.
    mats = struct('A', A, 'B', B, 'Q', Q, 'M', M, 'R', R);
else
    % Without the objective, things are much simpler.
    [i, C] = getindices([Nx, Nu]);
    C(i{1},i{1}) = args.a;
    C(i{1},i{2}) = args.b;
    M = expm(args.Delta*C);
    mats.A = M(i{1},i{1});
    mats.B = M(i{1},i{2});
    
    % Do quadrature for objective terms if necessary.
    if args.quad
        D = [args.q, args.m; args.m', args.r];
        integral = zeros(Nx + Nu, Nx + Nu);
        dt = args.Delta/(2*args.Nquad);
        for j = 0:2*args.Nquad
            if j == 0 || j == 2*args.Nquad
                mult = 1;
            else
                mult = 2 + 2*mod(j, 2);
            end
            expC = expm(j*dt*C);
            integral = integral + mult*expC'*D*expC;
        end
        integral = args.Delta/(6*args.Nquad)*integral;
        mats.Q = integral(i{1},i{1});
        mats.M = integral(i{1},i{2});
        mats.R = integral(i{2},i{2});
    end
end

% Handle linear terms.
if ~isempty(args.g) || ~isempty(args.h)
    if isempty(args.g)
        args.g = zeros(Nx, 1);
    end
    if isempty(args.h)
        args.h = zeros(Nu, 1);
    end
    [i, E] = getindices([Nx, Nu, Nx]);
    E(i{1},i{1}) = args.a;
    E(i{1},i{2}) = args.b;
    E(i{3},i{1}) = eye(Nx);
    E = expm(args.Delta*E);
    mats.G = E(i{3},i{1})'*args.g;
    mats.H = E(i{3},i{2})'*args.g + args.Delta*args.h;
end

% Choose return value.
if isequal(args.return, 'struct')
    varargout = {mats};
else
    varargout = cell(length(args.return), 1);
    fields = num2cell(args.return);
    [varargout{:}] = mpctools.structdeal(mats, fields{:});
end

end%function

function [r, nsquare] = expmnosquare(A)
    % [r, nsquare] = expmnosquare(A)
    %
    % Returns r and nsquare such that expm(A) = r^(2^nsquare). This is useful
    % if you have recursive formulas to preform the squaring yourself, e.g.,
    % as in Computing Integrals Involving the Matrix Exponential"
    % (Van Loan, 1987).
    %
    % Base implementation is taken from Octave but without trace rescaling or
    % balancing.
    narginchk(1, 1);
    n = size(A, 1);
    
    % Scaling.
    [~, e] = log2(norm(A, 'inf'));
    nsquare = min(max(0, e), 1023);
    a = A*2^(-nsquare);

    % Pade approximation for exp(A).
    c = [5.0000000000000000e-1,...
         1.1666666666666667e-1,...
         1.6666666666666667e-2,...
         1.6025641025641026e-3,...
         1.0683760683760684e-4,...
         4.8562548562548563e-6,...
         1.3875013875013875e-7,...
         1.9270852604185938e-9];

    a2 = a^2;
    id = eye(n);
    x = (((c(8)*a2 + c(6)*id)*a2 + c(4)*id)*a2 + c(2)*id)*a2 + id;
    y = (((c(7)*a2 + c(5)*id)*a2 + c(3)*id)*a2 + c(1)*id)*a;

    r = (x - y)\(x + y);
end%function

function [i, M] = getindices(blocks)
    % Returns block indices as a cell array i. Also returns zero matrix M.
    i = cell(length(blocks), 1);
    for j = 1:length(i)
        i0 = sum(blocks(1:(j - 1)));
        i{j} = (i0 + 1):(i0 + blocks(j));
    end
    M = zeros(sum(blocks));
end%function

