function num = varstruct2numstruct(var, fillval, varndims)
% num = varstruct2numstruct(var, fillval, [varndims])
%
% Creates a numeric struct num that is the same size as the variable struct var.
%
% varndims is a struct that gives the number of dimensions for each cell array.
% For any fields not given, the value is inferred from its size.
narginchk(2, 3);
if ~isscalar(fillval)
    error('fillval must be a scalar!');
end
if nargin() < 3
    varndims = struct();
end

% Loop through variable struct and add fields.
num = struct();
varnames = fieldnames(var);
for i = 1:length(varnames)
    v = varnames{i};
    
    % Short-circuit if variable is time-invariant.
    if ~iscell(var.(v))
        num.(v) = repmat(fillval, size(var.(v)));
        continue
    end
    
    % Check number of dimensions.
    cellsize = size(var.(v));
    nd = mpctools.structget(varndims, v, []);
    if isempty(nd)
        if cellsize(1) == 1
            nd = 1;
        else
            nd = length(cellsize);
        end
    end
    
    % Pad or clip as necessary to get correct number of dimensions.
    if nd > length(cellsize)
        cellsize = [cellsize, ones(1, nd - length(cellsize))];
    elseif nd < length(cellsize)
        front = 1:(length(cellsize) - nd);
        nbad = sum(cellsize(front) ~= 1);
        if nbad > 0
            error('Field "%s" is supposed to have %d dimensions but has %d!', ...
                  v, nd, nd + nbad);
        end
        cellsize(front) = [];
    end
    
    % Fill with fill value.
    vecsize = var.(v){1}.size;
    if vecsize(2) ~= 1
        error('Variable "%s" is not a vector!', v);
    end
    thissize = [vecsize(1), cellsize];
    num.(v) = repmat(fillval, thissize);
end
end%function

