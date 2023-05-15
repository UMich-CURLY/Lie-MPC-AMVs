function abscon = absconstraints(var, name)
% abscon = absconstraints(var, name)
%
% Returns a list of constraints to define absolute values, i.e., with name 'x',
%
%    absx >= x,   absx >= -x
%
% with one constraint written for each element of x.
if ~all(isfield(var, {name, ['abs', name]}))
    error('var must have fields "%s" and "abs%s"!', name, name);
end

x = var.(name);
absx = var.(['abs', name]);
if ~isequal(size(x), size(absx))
    error('Variables must be equal sizes!');
end

consize = size(x);
con = cell([2, consize]);
for i = 1:prod(consize)
    con{2*i - 1} = x{i} - absx{i};
    con{2*i} = -x{i} - absx{i};
end
Nx = size(x{1}, 1);
abscon = struct();
abscon.con = con;
abscon.lb = -inf([Nx, 2, consize]);
abscon.ub = zeros(size(abscon.lb));

end%function

