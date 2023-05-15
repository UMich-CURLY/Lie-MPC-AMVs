function [shapes, repeats] = defaultvarsizes_(Nt, finalx, finaly)
% [shapes, repeats] = defaultvarsizes_(Nt, finalx, finaly)
%
% Returns a struct of variable sizes and repeats.
% The second and third arguments are optional and default to True.
narginchk(0, 3);

% Default values.
if nargin() < 1 || isempty(Nt)
    Nt = 1;
end
if nargin() < 2 || isempty(finalx)
    finalx = true();
end
if nargin() < 3 || isempty(finaly)
    finaly = true();
end

% Create a struct of all possible variables. First entry is vector shape, second
% entry is number of repeats.
sizes = struct();
sizes.x = {'x', Nt + finalx};
sizes.z = {'z', Nt + finalx};
sizes.u = {'u', Nt};
sizes.y = {'y', Nt + finaly};
sizes.w = {'w', Nt};
sizes.v = {'v', Nt + finaly};
sizes.xc = {'x', {'c', Nt}};
sizes.zc = {'z', {'c', Nt}};
sizes.s = {'s', Nt};
sizes.d = {'d', Nt + (finalx || finaly)};

sizes.xsp = sizes.x;
sizes.usp = sizes.u;
sizes.Dx = {'x', Nt - 1 + finalx};
sizes.Du = sizes.u;
sizes.Dd = {'d', Nt - 1 + (finalx || finaly)};

absvars = {'x', 'u', 'Dx', 'Du', 'xc', 'w', 'v', 'y', 'd'};
for i = 1:length(absvars)
    v = absvars{i};
    absv = ['abs', v];
    sizes.(absv) = sizes.(v);
end

% Now split into separate shapes and repeats.
shapes = struct();
repeats = struct();
varnames = fieldnames(sizes);
for i = 1:length(varnames)
    v = varnames{i};
    s = sizes.(v);
    shapes.(v) = promoteshape(s{1});
    repeats.(v) = promoterepeat(s{2});
end

end%function

function s = promoteshape(s)
    % Promotes scalar shape to cell array as necessary.
    if ~iscell(s)
        s = {s, 1};
    end
end%function

function s = promoterepeat(s)
    % Promotes a scalar repeat to a cell array as necessary. Note that no
    % padding is added, in contrast to promoteshape.
    if ~iscell(s)
        s = {s};
    end
end%function

