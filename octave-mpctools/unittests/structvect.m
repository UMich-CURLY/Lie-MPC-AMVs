function v = structvect(s)
% v = structvect(s)
%
% Concatenates the fields of the struct s into a single vector, accounting for
% time.
%
% Assumes all fields have time running along the final dimension. This ensures
% that variables are concatenated as `[x_0; u_0; x_1; u_1]`, etc.
%
% Note: this is a reference implementation for the "simple" case where
% everything is assumed to be time-varying. It is used in the unit test to
% compare the new implementation, but it is not used inside actual mpctools
% functions.
fields = mpctools.VarLayout.sortnames(fieldnames(s));
Nf = length(fields);
Nts = structfun(@(f) size(f, ndims(f)), s, 'UniformOutput', false());
Nt = max(structfun(@(f) size(f, ndims(f)), s));
vecs = cell(Nf, Nt);
keep = true(size(vecs));
for i = 1:Nf
    f = fields{i};
    for t = 1:Nts.(f)
        val = indexfinaldim(s.(f), t);
        if iscell(val)
            val = flattencell(val);
        else
            val = val(:);
        end
        vecs{i, t} = val;
    end
    keep(i,Nts.(f) + 1:end) = false();
end
vecs = vecs(keep(:));
v = vertcat(vecs{:});

end%function

function xt = indexfinaldim(x, t)
    % `xt = indexfinaldim(x, t)`
    %
    % Returns `x(...,t)`, where `...` is an appropriate number of colons.
    switch ndims(x)
    case 2
        xt = x(:,t);
    case 3
        xt = x(:,:,t);
    otherwise
        idx = substruct('()', [repmat({':'}, 1, ndims(x) - 1), {t}]);
        xt = subsref(x, idx);     
    end
end%function

function cflat = flattencell(c)
    % Flattens each element of the cell array c and then concatenates.
    for i = 1:numel(c)
        c{i} = c{i}(:);
    end
    cflat = vertcat(c{:});
end%function

