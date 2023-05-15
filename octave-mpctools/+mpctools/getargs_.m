function args = getargs_(idx, vars, names)
% args = getargs_(idx, vars, names)
%
% Returns cell array of arguments from fields of struct vars in the order
% listed in names.
%
% Each element of vars should be a cell array of variables that is indexed
% using the rightmost entries of the index vector i. Note that this means any
% cell arrays that are supposed to have only one dimension must be row vectors.
%
% All cell arrays are accessed modulo length (along each dimension). For
% example, for a cell array of length 2, if idx is 0, the second entry is
% returned, while is idx is 3, the first entry is returned.
%
% Note that names can also be a casadi.Function, at which point the names are
% read from the function's inputs.
narginchk(3, 3);
if ~iscell(names)
    if mpctools.iscasadifunc(names)
        names = mpctools.names2cell_(names.name_in());
    else
        error('Invalid input for names!');
    end
end
args = cell(1, length(names));
for i = 1:length(args)
    n = names{i};
    v = vars.(n);
    
    if ~iscell(v)
        % Shortcut if variable is "scalar" (i.e., only has one copy).
        args{i} = v;
    else
        % Check number of dimensions.
        if isrow(v)
            ni = 1;
        else
            ni = sum(size(v) > 1);
        end
        if ni > length(idx)
            error('Entry "%s" takes %d indices, but only %d provided!', ...
                  n, ni, length(idx))
        end
        
        % Do the actual indexing.
        if ni == 1
            args{i} = cellmod1d(v, idx(end));
        else
            if length(idx) < ni
                error('Variable "%s" takes %d indices but only %d given!', ...
                      n, ni, length(idx));
            end
            args{i} = cellmodNd(vars.(n), idx((end - ni + 1):end));
        end
    end
end

end%function

function ci = cellmod1d(c, i)
    % Returns c{i mod length(c)}
    i = mod(i - 1, length(c)) + 1;
    ci = c{i};
end%function

function ci = cellmodNd(c, i)
    % Returns c{i mod size(c)} with i multidimensional.
    shape = size(c);
    i = mod(i - 1, shape) + 1;
    ci = subsref(c, substruct('{}', num2cell(i)));
end%function

