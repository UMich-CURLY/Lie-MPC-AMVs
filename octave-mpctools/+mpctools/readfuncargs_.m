function funcargs = readfuncargs_(funcnames, funcs)
% funcargs = readfuncargs_(funcnames, funcs)
%
% For each function name "f" listed in funcnames, checks if funcs.f is given,
% and if so, adds its arguments to funcargs.f.
%
% funcs should be a struct of Casadi Function objects.
narginchk(2, 2);
funcargs = struct();
for i = 1:length(funcnames)
    f = funcnames{i};
    if isfield(funcs, f) && ~isempty(funcs.(f))
        funcargs.(f) = mpctools.names2cell_(funcs.(f).name_in());
    end
end

end%function

