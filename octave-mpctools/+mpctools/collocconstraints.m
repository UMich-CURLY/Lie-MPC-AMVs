function model = collocconstraints(f, var, funcargs, A, Delta)
% model = collocconstraints(f, var, funcargs, A, Delta)
%
% Returns a struct with field 'con' containing cell array of constraints for
% orthogonal collocation as well as fields 'lb' and 'ub' with vectors of bounds.

% Get sizes.
if ~isfield(var, 'x')
    error('var must have an "x" entry!');
end
Nc = size(var.x, 1) - 2;
Nt = size(var.x, 2);
Nx = size(var.x{1}, 1);

% Make the constraints.
cons = cell(Nc + 1, Nt);
for k = 1:Nt
    X = horzcat(var.x{:,k});
    for i = 1:(Nc + 1)
        fargs = mpctools.getargs_([i, k], var, funcargs);
        try
            dxdt = Delta*f(fargs{:});
        catch err
            error('evaluating f(%s): %s', mpctools.row2str(funcargs), ...
                  err.message);
        end
        cons{i,k} = dxdt - X*A(i,:)';
    end
end

% Assemble into struct.
model = struct();
model.con = cons;
model.lb = zeros(Nx, Nc + 1, Nt);
model.ub = model.lb;
end%function

