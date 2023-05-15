function varandpar = mergevarandpar_(var, par, allargs)
% varandpar = mergevarandpar_(var, par, allargs)
%
% Merges variables and parameters into a single struct and makes sure that each
% element of allargs is represented.

% Check for overlaps in parameters and variables.
clashes = intersect(fieldnames(var), fieldnames(par));
if ~isempty(clashes)
    fmt = repmat(' "%s"', 1, length(clashes));
    error(['Clashing parameter names:', fmt], clashes{:});
end
varandpar = mpctools.structupdate(var, par);

% Make sure each element of allargs is in varandpar.
unknownargs = setdiff(allargs, fieldnames(varandpar));
if ~isempty(unknownargs)
    if ismember('s', unknownargs)
        hint = ' (did you forget to specify N.s?)';
    else
        hint = '';
    end
    fmt = repmat(' "%s"', 1, length(unknownargs));
    error(['Unknown function arguments:', fmt, hint], unknownargs{:});
end

end%function
