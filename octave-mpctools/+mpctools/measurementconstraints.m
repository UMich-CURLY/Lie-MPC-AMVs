function model = measurementconstraints(h, vars, funcargs)
% model = measurementconstraints(h, vars, funcargs)
%
% Returns a struct with field 'con' containing a cell array of measurement
% constraints and fields 'lb' and 'ub' containing lower and upper bound vectors.
if ~isfield(vars, 'y')
    error('Field "y" is missing from vars!');
end
y = vars.y;
Nt = length(y);
Ny = size(y{1}, 1);
cons = cell(1, Nt);
for i = 1:Nt
    args = mpctools.getargs_(i, vars, funcargs);
    try
        cons{i} = h(args{:}) - y{i};
    catch err
        error('evaluating h(%s): %s', mpctools.row2str(funcargs), err.message);
    end
end

model = struct();
model.con = cons;
model.lb = zeros(Ny, Nt);
model.ub = model.lb;
end%function

