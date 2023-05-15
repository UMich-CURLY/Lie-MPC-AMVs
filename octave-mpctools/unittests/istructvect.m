function s = istructvect(v, sizes)
% s = istructvect(v, sizes)
%
% Splits the large vector v into appropriate struct fields, respecting time.
% Inverse operation of structvect.
%
% Assumes all fields have time running along the final dimension, i.e., that
% variables were concatenated as `[x_0; u_0; x_1; u_1]` etc.
%
% Note: this is a reference implementation for the "simple" case where
% everything is assumed to be time-varying. It is used in the unit test to
% compare the new implementation, but it is not used inside actual mpctools
% functions.
if ~isnumeric(v)
    error('Input must be numeric!');
end
fields = mpctools.VarLayout.sortnames(fieldnames(sizes));
Nf = length(fields);
v = v(:);
Nts = structfun(@(s) s(end), sizes, 'UniformOutput', false());
Nt = max(structfun(@(s) s(end), sizes));

% Preallocate cell arrays.
s = struct();
for i = 1:Nf
    f = fields{i};
    s.(f) = cell(Nts.(f), 1);
end

% Loop through time and fields.
row = 0;
for t = 1:Nt
    for i = 1:Nf
        f = fields{i};
        if t <= Nts.(f)
            chunksize = sizes.(f)(1:end-1);
            if isscalar(chunksize)
                chunksize = [chunksize, 1];
            end
            chunk = reshape(v((row + 1):(row + prod(chunksize))), chunksize);
            s.(f){t} = chunk;
            row = row + prod(chunksize);
        end
    end
end

% Concatenate chunks.
for i = 1:Nf
    f = fields{i};
    s.(f) = cat(length(sizes.(f)), s.(f){:});
end

end%function

