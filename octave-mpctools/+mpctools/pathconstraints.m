function pathcon = pathconstraints(Nt, e, vars, funcargs)
% pathcon = pathconstraints(Nt, e, vars, funcargs)
%
% Formulates nonlinear path constraints.
cons = cell(1, Nt);
for i = 1:Nt
    args = mpctools.getargs_(i, vars, funcargs);
    try
        cons{i} = e(args{:});
    catch err
        error('evaluating e(%s): %s', mpctools.row2str(funcargs), err.message);
    end
    if i == 1
        Ne = size(cons{i}, 1);
    end
end

% Package into struct.
pathcon = struct();
pathcon.con = cons;
pathcon.lb = -inf(Ne, Nt);
pathcon.ub = zeros(Ne, Nt);

end%function

