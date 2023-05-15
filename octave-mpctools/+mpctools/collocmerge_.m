function vars = collocmerge_(vars, funcargs, jumpvars)
% var = collocmerge_(vars, funcargs, [jumpvars])
%
% Merges the 'x' and 'xc' fields of vars. Also merges 'z', and 'zc' if present.
%
% If jumpvars is given, it should be a struct with fields 'x' and/or 'u' giving
% variables that give the "jump" between time points (e.g., w in MHE). The jumps
% are assumed to happen at the end of each collocation interval.
narginchk(2, 3);
if nargin() < 3
    jumpvars = struct();
end
for v = ['x', 'z']
    if isfield(vars, v)
        vars = collocmerge(vars, jumpvars, funcargs, v);
    end
end

end%function

function vars = collocmerge(vars, jumpvars, funcargs, name)
    % name should be a string 'x' or 'z' to indicate which variable is being
    % used.
    namec = [name, 'c'];
    
    if ~all(isfield(vars, {name, namec}))
        error('Variables %s and %s are required for collocation!', ...
              name, namec);
    elseif any(structfun(@(a) ismember(namec, a), funcargs))
        error('%s cannot appear in funcargs (replace with %s)!', ...
              namec, name);
    end
    
    v = vars.(name);
    vc = vars.(namec);
    vjump = mpctools.structget(jumpvars, name, []);
    Nt = length(v) - 1;
    if size(vc, 2) ~= Nt
        error('%s has the wrong number of time points!', namec);
    end
    Nc = size(vc, 1);
    Nv = size(v{1}, 1);
    if ~isempty(vjump)
        if length(vjump) ~= Nt
            error('%s jump variable has the wrong number of time points!', name);
        elseif size(vjump{1}, 1) ~= Nv
            error('%s jump variable has the wrong number of components!', name);
        end
    end

    % Duplicate the endpoints in v and absv.
    vars = rmfield(vars, {name, namec});
    vars.(name) = dupep(v, vc, vjump, Nc, Nt);
    absvars = {['abs', name], ['abs', namec]};
    if all(isfield(vars, absvars))
        absv = vars.(absvars{1});
        absvc = vars.(absvars{2});
        vars = rmfield(vars, absvars);
        vars.(absvars{1}) = dupep(absv, absvc, Nc, Nt);
    end
end%function

function X = dupep(x, xc, xjump, Nc, Nt)
    % Duplicates the x variables as endpoints of xc.
    X = cell(Nc + 2, Nt);
    for k = 1:Nt
        X(:,k) = [x(k); xc(:,k); x(k + 1)];
        if ~isempty(xjump)
            X{end,k} = X{end,k} - xjump{k};
        end
    end
end%function

