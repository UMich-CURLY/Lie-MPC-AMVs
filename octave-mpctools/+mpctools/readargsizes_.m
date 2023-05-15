function [argsizes, allargs] = readargsizes_(funcargs, funcs)
% [argsizes, allargs] = readargsizes_(funcargs, funcs)
%
% Returns sizes of all arguments defined in funcargs.
%
% funcargs should be a struct of function names and variable names. funcs should
% be a struct of casadi.Function objects.
%
% argsizes is a struct of argument sizes. Note that a given argument must be the
% same size in all functions. allargs is a list of all the arguments that were
% found (as a row cell array).
narginchk(2, 2);
funcnames = fieldnames(funcargs);
allargs = {};
argsizes = struct();
wherefound = struct();
for i = 1:length(funcnames)
    f = funcnames{i};
    func = funcs.(f);
    fargs = funcargs.(f);
    for j = 1:length(fargs)
        name = fargs{j};
        shape = func.size_in(j - 1); % Note 0-based indexing here.
        if isfield(argsizes, name)
            checkshape = argsizes.(name);
            if ~isequal(shape, checkshape)
                s1 = shapemessage(f, shape);
                s2 = shapemessage(wherefound.(name), checkshape);
                error('Argument %s is %s but %s!', name, s1, s2);
            end
        else
            argsizes.(name) = shape;
            wherefound.(name) = f;
        end
    end
end
allargs = fieldnames(argsizes)';

end%function

function m = shapemessage(f, s)
    % Returns a message string for the given function and shape.
    m = sprintf('shape [%s] in function %s', mpctools.row2str(s), f);
end

