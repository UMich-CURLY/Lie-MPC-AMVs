function [algebra, collocalg] = algebraconstraints(g, vars, funcargs, varsc)
% [algebra, collocalg] = algebraconstraints(g, vars, funcargs, [varsc])
%
% Returns a struct with field 'con' containing a cell array of algebra
% constraints and fields 'lb' and 'ub' containing lower and upper bound vectors.
%
% Note that if g takes 'u' as an argument and `vars.u` has only `N - 1`
% elements, then `u(N - 1)` is used in the final algebraic constraint.
% That is, for time `t = N`, the constraint becomes
%
%     $g(x(N), z(N), u(N - 1)) = 0$
%
% If you do not want this behavior, you should add the appropriate final
% element so that `vars.u` has `N` elements like `vars.x` and `vars.u`.
%
% If varsc is provided, constraints are also returned for the interior
% collocation points.
narginchk(3, 4);

% Get sizes.
if ~isfield(vars, 'z')
    error('vars must have a "z" entry!');
end

% Special case for u variables.
Nt = size(vars.z, 2);
if ismember('u', funcargs) && isfield(vars, 'u') && length(vars.u) == Nt - 1
    vars.u = [vars.u, vars.u(end)];
end

% Constraints for boundary points.
algebra = struct();
algebra.con = docon(g, vars, funcargs, 1, 1:size(vars.z, 2));

Ng = size(algebra.con{1,1}, 1);
algebra.lb = zeros(Ng, Nt);
algebra.ub = algebra.lb;

% Constraints for interior collocation points.
if nargin() >= 4
    if ~isfield(varsc, 'z')
        error('varsc must have a "z" entry!');
    end
    collocalg = struct();
    [Nc, Nt] = size(varsc.z);
    collocalg.con = docon(g, varsc, funcargs, 2:(Nc - 1), 1:Nt);
    collocalg.lb = zeros(Ng, Nc - 2, Nt);
    collocalg.ub = collocalg.lb;
else
    collocalg = mpctools.emptycon();
end

end%function

function cons = docon(g, vars, funcargs, Cinds, Tinds)
    cons = cell(Cinds(end), Tinds(end));
    for k = Tinds
        for i = Cinds
            gargs = mpctools.getargs_([i, k], vars, funcargs);
            try
                cons{i,k} = g(gargs{:});
            catch err
                error('evaluating g(%s): %s', mpctools.row2str(funcargs), ...
                      err.message);
            end
        end
    end
    cons = cons(Cinds, Tinds);
end%function

