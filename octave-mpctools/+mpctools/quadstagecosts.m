function obj = quadstagecosts(Nt, l, var, funcargs, obj, q, Delta)
% obj = quadstagecosts(Nt, l, var, funcargs, obj, q, Delta)
%
% Returns an SX object for the objective function via collocation-based
% quadrature.
narginchk(7, 7);
Nc = length(q);
terms = cell(Nc, Nt);
for t = 1:Nt
    for c = 1:Nc
        args = mpctools.getargs_([c, t], var, funcargs);
        try
            terms{c, t} = q(c)*Delta*l(args{:});
        catch err
            error('evaluating l(%s): %s', mpctools.row2str(funcargs), ...
                  err.message);
        end
    end
end
terms = vertcat(terms{:,:});
obj = obj + terms.sum();
end%function

