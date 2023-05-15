function varargout = collocweights(N, varargin)
% `[r, a, b, q] = collocweights(N, ['left'], ['right'])`
%
% Returns collocation weights based on roots of Legendre polynomials.
%
% The first argument `N` gives the number of interior collocation points. The
% second and third arguments can be the strings 'left' or 'right' to add
% additional collocation points on the boundary.
%
% Note that all arguments must be passed as positional arguments.
%
% In Octave, simply calls builtin colloc. In Matlab, uses m-file implementation.
%
% Reference: J. Villadsen, M. L. Michelsen, "Solution of Differential
%            Equation Models by Polynomial Approximation".
varargout = cell(1, nargout());
if mpctools.isOctave() && ~ismember('matlab', varargin)
    [varargout{:}] = colloc(N, varargin{:});
else
    [varargout{:}] = colloc_matlab(N, varargin{:});
end

end%function

% *****************************************************************************
% m-file implementation for Matlab below.
% *****************************************************************************

function [r, a, b, q] = colloc_matlab(N, varargin)
% m-file implementation of Octave's colloc() for Matlab.
narginchk(1, 3);
if ~isscalar(N) || N < 0 || round(N) ~= N
    error('N must be a nonnegative integer!');
end

% Use Casadi's collocation_points() function to get roots if N < 10. Otherwise,
% use a matrix method.
if N < 10
    r = casadi.collocation_points(N, 'legendre');
else
    r = legendreroots(N);
end

% Add endpoint if requested.
if ismember('left', varargin)
    r = [0, r];
    N = N + 1;
end
if ismember('right', varargin)
    r = [r, 1];
    N = N + 1;
end

% Calculate other guys.
if nargout() > 1
    % Calculate derivatives.
    d1 = ones(1, N);
    d2 = zeros(1, N);
    d3 = zeros(1, N);

    % Use recursive formulas.
    for j = 1:N
        dr = r - r(j);
        dr(j) = 1;
        mul2 = 3*ones(size(d2));
        mul2(j) = 0;
        mul1 = 2*ones(size(d1));
        mul1(j) = 0;

        d3 = dr.*d3 + mul2.*d2;
        d2 = dr.*d2 + mul1.*d1;
        d1 = dr.*d1;
    end
    
    a = dfopr(N, d1, d2, d3, r, 'first');
    if nargout() > 2
        b = dfopr(N, d1, d2, d3, r, 'second');
        if nargout() > 3
            q = dfopr(N, d1, d2, d3, r, 'weights');
            q = q(:); % Column vector.
        end
    end
end

r = r(:); % Column vector.

end%function

function z = bincoef(n, k)
    % Returns binomial coefficient for n and k with extrapolation to non-
    % integral values.
    %
    % Not calculated very carefully, so don't use this for large values of
    % n or k.
    z = gamma(n + 1)./gamma(k + 1)./gamma(n - k + 1);
end%function

function M = dfopr(N, d1, d2, d3, r, matrix)
    % Calculates "weights", "first", or "second" derivative matrices.
    narginchk(6, 6);
    switch matrix
    case 'weights'
        ax = r.*(1 - r);
        
        % Check if 0 and 1 are endpoints.
        if r(1) ~= 0
            ax = ax./r.^2;
        end
        if r(end) ~= 1
            ax = ax./(1 - r).^2;
        end
        M = ax./d1.^2;
        M = M/sum(M(:));
    case {'first', 'second'}
        % First handle diagonal elements.
        if isequal(matrix, 'first')
            m = d2./(2*d1);
        else
            m = d3./(3*d1);
        end
        
        % Now do off-diagonal stuff.
        [ri, rj] = ndgrid(r, r);
        [d1i, d1j] = ndgrid(d1, d1);
        [d2i, d2j] = ndgrid(d2, d2);
        
        y = ri - rj;
        diagInds = (0:(N - 1))*N + (1:N); % Diagnoal indices.
        y(diagInds) = 1; % We will fix this later.
        M = d1i./(d1j.*y);
        if isequal(matrix, 'second')
            M = M.*(d2i./d1i - 2./y);
        end
        M(diagInds) = m; % Add back in diagonal elements.
    otherwise
        error('Unknown value: "%s"', matrix);
    end
end%function

function r = legendreroots(N)
    % r = legendreroots(N)
    %
    % Returns the N distinct roots of the Nth-order Legendre polynomial by
    % finding the eigenvalues of the Jacobi matrix associated with the
    % recurrence relationship
    %
    %    P_{n + 1}(x) = x P_n(x) - n^2/(4*n^2 - 1) P_{n - 1}(x)
    %
    % Note that the matrix is shifted so that the roots lie in (0, 1) rather
    % than (-1, 1).
    %
    % This function may be inaccurate for extremely large values of 
    n = (1:(N - 1))';
    gam = n./sqrt(4*n.^2 - 1);
    A = 0.5*spdiags([[0; gam], [gam; 0]], [1, -1], speye(N));
    r = eig(A)';
end%function

