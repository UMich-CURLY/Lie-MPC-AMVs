function roc = rocconstraints(var, whichvar)
% roc = rocconstraints(var, 'u')
% roc = rocconstraints(var, 'x')
% roc = rocconstraints(var, 'd')
%
% Returns a struct with field 'con' containing a cell array of rate-of-change
% constraints and fields 'lb' and 'ub' containing lower and upper bound vectors.
%
% For 'u' constraints, var must contain 'u', 'Du', and 'uprev'. For 'x',
% constraints, only 'x', and 'Dx' are required.
narginchk(2, 2);
switch whichvar
case 'u'
    if ~all(isfield(var, {'u', 'Du', 'uprev'}))
        error('var must contain fields "u", "Du", and "uprev"!');
    end
    u = var.u;
    Du = var.Du;
    uprev = var.uprev;
case 'x'
    if ~all(isfield(var, {'x', 'Dx'}))
        error('var must contain fields "x" and "Dx"!');
    end
    u = var.x(2:end);
    Du = var.Dx;
    uprev = var.x{1};
case 'd'
    if ~all(isfield(var, {'d', 'Dd'}))
        error('var must contain fields "d" and "Dd"!');
    end
    u = var.d(2:end);
    Du = var.Dd;
    uprev = var.d{1};
otherwise
    error('Invalid variable: %s', whichvar);
end
% TODO: remove switch and just check for whichvar, and ['D', whichvar] fields.

% Internally, we call everything u, even though it could be x or u.
Nu = size(uprev, 1);
Nt = length(u);
if length(Du) ~= Nt
    error('Incorrect length for var.D%s!', whichvar);
end
cons = cell(1, Nt);
cons{1} = u{1} - uprev - Du{1};
for i = 2:Nt
    cons{i} = u{i} - u{i - 1} - Du{i};
end
roc = struct();
roc.con = cons;
roc.lb = zeros(Nu, Nt);
roc.ub = roc.lb;
end%function

