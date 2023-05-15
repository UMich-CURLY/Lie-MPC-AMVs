function [I, name] = getargindex_(func, K, inorout)
% [i, name] = getargindex_(func, k)
% [i, name] = getargindex_(func, k, 'in')
% [i, name] = getargindex_(func, k, 'out')
% [i, name] = getargindex_(argcell, k, 'cell')
%
% Returns the (1-based) integer index and argument name of CasADi function
% `func`. `k` can be either the integer index, or a string giving the variable
% name.
%
% The third argument specifies whether `k` is a function input or output.
%
% `k` can also be a cell array, at which point i will be an array, and name
% will be a cell array with the same size as k.
%
% In the final syntax, the first argument should be a row cell array of
% argument names.
narginchk(2, 3);
if nargin() < 3
    inorout = 'in';
end

% Check for cell array input.
kcell = iscell(K);
if ~kcell
    K = {K};
end

% Get names of arguments.
switch inorout
case {'in', 'input'}
    names = mpctools.names2cell_(func.name_in());
    inorout = 'input';
    funcname = ['function ', func.name()];
case {'out', 'output'}
    names = mpctools.names2cell_(func.name_out());
    inorout = 'output';
    funcname = ['function ', func.name()];
case 'cell'
    if ~iscell(func)
        error('First argument must be a cell array to use ''cell''!');
    end
    names = func(:)';
    inorout = 'argument';
    funcname = mpctools.row2str(names);
otherwise
    error('Invalid third argument!');
end

% Loop through elements of K and decide for each.
I = NaN(size(K));
name = cell(size(K));
for j = 1:numel(K)
    k = K{j};
    if ischar(k)
        [~, i] = ismember(k, names);
        if isempty(i)
            error('Name %s is not an %s of %s!', k, inorout, funcname);
        end
        I(j) = i;
        name{j} = k;
    elseif isscalar(k) && round(k) == k && k > 0
        i = k;
        if i > length(names)
            error('Index %d is too large: %s has only %d %ss!', ...
                  i, funcname, length(names), inorout);
        end
        I(j) = i;
        name{j} = names{i};
    else
        error('Unknown input for k! Must be an integer or string!');
    end
end

% Remove cell array from name if necessary.
if ~kcell
    name = name{1};
end

end%function

