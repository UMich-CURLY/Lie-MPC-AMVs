function var = getSymStruct(casaditype, sizes, repeats)
% var = getSymStruct(casaditype, sizes, repeats)
%
% Returns a struct of cell arrays of Casadi SX or MX variables.
%
% sizes should be a struct of sizes for the individual symbolic variables.
% repeats should be a struct giving the number of copies (default is a single
% copy that isn't in a cell array).

% Check arguments.
narginchk(3, 3);
if ~isstruct(sizes) || ~isstruct(repeats)
    error('sizes and repeats must be structs!');
end

% Get names and convert other arguments to cell arrays.
names = fieldnames(sizes);
sizes = cellfun(@(f) sizes.(f), names, 'UniformOutput', false());
repeats = cellfun(@(f) mpctools.structget(repeats, f, []), names, ...
                  'UniformOutput', false());

% Make struct.
var = getstruct_(names, sizes, repeats, casaditype);

end%function

function vars = getstruct_(names, vecsizes, cellsizes, casaditype)
    % Assembles the struct using casadi symbolic (SX or MX) variables.
    narginchk(4, 4);
    symfunc = mpctools.getcasadisymfunc(casaditype);
    vars = struct();
    for n = 1:length(names)
        name = names{n};
        thisvecsize = vecsizes{n};
        strcellsize = cellsizes{n};
        thiscellsize = strcellsize;
        if isscalar(thiscellsize)
            thiscellsize = [1, thiscellsize];
        end
        if isempty(thiscellsize)
            % Don't need cell array.
            vars.(name) = symfunc(name, thisvecsize);
        else
            % Do need cell array.
            pad = ones(1, 2 - length(thiscellsize));
            thiscell = cell([pad, thiscellsize]);
            for i = 1:prod(thiscellsize)
                thiscell{i} = symvec(i, name, thisvecsize, strcellsize, symfunc);
            end
            vars.(name) = thiscell;
        end
    end
end%function

function x = symvec(i, name, vecsize, cellsize, symfunc)
    % Returns symbolic vector of appropriate size.
    name = sprintf('%s%s', name, ind2str(i, cellsize));
    x = symfunc(name, vecsize);
end%function

function s = ind2str(ind, varsize)
    % s = ind2str(ind, varsize)
    %
    % Returns string representation of subscripts.
    if length(varsize) == 1
        s = sprintf('(%d)', ind - 1); % 0-based index.
    else
        subs = cell(length(varsize), 1);
        [subs{:}] = ind2sub(varsize, ind);
        subs = cell2mat(subs) - 1; % Convert to 0-based index.
        s = ['(', mpctools.row2str(subs), ')'];
    end
end%function

