function [X, T, xcflat, tflat] = smushcolloc(x, xc, t, tc)
% [X, T, xcflat, tcflat] = smushcolloc(x, xc, [t], [tc])
%
% Combines [Nx by Nt + 1] x and [Nx by Nc by Nt] xc into a single
% [Nx by (Nc + 1)*Nt + 1] array.
%
% Also combines t and tc if given (otherwise, they are calculated). Note that
% t must be a vector of length Nt + 1, and tc must be a [Nc by Nt] matrix.
narginchk(2, 4);

% Check x and xc.
Nx = size(xc, 1);
Nc = size(xc, 2);
Nt = size(xc, 3);
if ~isequal(size(x), [Nx, Nt + 1])
    error('Incorrect size for x!');
end

% Make X.
xc = reshape(xc, [Nx*Nc, Nt]);
X = [x(:,1:end-1); xc];
X = [reshape(X, [Nx, (Nc + 1)*Nt]), x(:,end)];

% Make T if needed.
if nargout() > 1
    % Check t.
    if nargin() < 3 || isempty(t)
        t = 0:Nt;
    elseif ~isvector(t) || length(t) ~= Nt + 1
        error('Incorrect size for t!');
    end
    t = t(:)'; % Make sure it's a row vector.

    % Check tc.
    if nargin() < 4 || isempty(tc)
        r = mpctools.collocweights(Nc);
        tc = bsxfun(@plus, t(1:end-1), bsxfun(@times, diff(t), r));
    elseif ~isequal(size(tc), [Nc, Nt])
        error('Incorrect size for tc!');
    end

    % Do the actual work.
    T = [t(1:end-1); tc];
    T = [T(:)', t(end)];
    
    % Flatten the other guys for easier plotting.
    xcflat = reshape(xc, Nx, Nc*Nt);
    tcflat = reshape(tc, 1, Nc*Nt);
end

end%function.
