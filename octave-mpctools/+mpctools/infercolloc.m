function xc = infercolloc(r, x, method, infval)
% xc = infercolloc(r, x, [method])
%
% Interpolates the `[Nx, Nt + 1]` array `x` at the collocation points defined in
% `r` to yield an `[Nx, Nc, Nt]` array of values at the collocation points.
%
% Default for `method` is 'linear'. Other options are 'pchip' and 'spline'. See
% `help interp1` for details.
narginchk(2, 4);
if nargin() < 3
    method = 'linear';
end
if nargin() < 4
    infval = 1e100;
end

Nx = size(x, 1);
Nt = size(x, 2) - 1;

r = r(2:end-1); % Get rid of endpoints.
r = r(:)'; % Make sure it's a row vector.
Nc = length(r);

if Nt == 0
    xc = repmat(x, [1, Nc, 1]);
else
    xinf = isinf(x);
    x(xinf) = infval*sign(x(xinf));
    
    t = 0:Nt;
    tc = kron(1:Nt, r);
    xc = interp1(t, x', tc)';

    xc = reshape(xc, [Nx, Nc, Nt]);
    xcinf = abs(xc) > sqrt(infval);
    xc(xcinf) = inf()*sign(xc(xcinf));
end

end%function

