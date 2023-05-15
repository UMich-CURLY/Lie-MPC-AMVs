function model = multipleshooting(f, vars, funcargs)
% model = multipleshooting(f, vars, funcargs)
%
% Returns a struct with field 'con' containing a cell array of multiple shooting
% constraints and fields 'lb' and 'ub' containing lower and upper bound vectors.
if ~isfield(vars, 'x')
    error('Field "x" is missing from vars!');
end
x = vars.x;
Nt = length(x) - 1;
Nx = size(x{1}, 1);
cons = cell(1, Nt);
for i = 1:Nt
    args = mpctools.getargs_(i, vars, funcargs);
    try
        cons{i} = f(args{:}) - x{i + 1};
    catch err
        error('evaluating f(%s): %s', mpctools.row2str(funcargs), err.message);
    end
end

model = struct();
model.con = cons;
model.lb = zeros(Nx, Nt);
model.ub = model.lb;
end%function

