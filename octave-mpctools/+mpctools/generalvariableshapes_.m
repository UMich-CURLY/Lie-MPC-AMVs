function [shapes, repeats, allshapes] = generalvariableshapes_(varlist, N, finalx, finaly)
% [varshapes, parshapes, allshapes] = generalvariableshapes_(varlist, N,
%                                                 [finalx=True], [finaly=True])
%
% Returns a struct of shapes and repeats for all variables.
narginchk(2, 4);
if nargin() < 3
    finalx = true();
end
if nargin() < 4
    finaly = true();
end

% Get number of time points.
if ~isfield(N, 't')
    error('N must contain a "t" entry!');
end
Nt = N.t;

% Get a list of all possible shapes.
[allshapes, allrepeats] = mpctools.defaultvarsizes_(Nt, finalx, finaly);

% Loop through variables.
shapes = struct();
repeats = struct();
for i = 1:length(varlist)
    v = varlist{i};
    if ~isfield(allshapes, v)
        error('Unknown variable name "%s"', v);
    end
    shapes.(v) = allshapes.(v);
    repeats.(v) = allrepeats.(v);
end

% Change cell arrays to row vectors.
shapes = cell2shape(shapes, N);
repeats = cell2shape(repeats, N);
allshapes = cell2shape(allshapes, N);

end%function

function out = cell2shape(shapes, N)
    % Converts cell array with ints and strings to row vector.
    out = struct();
    names = fieldnames(shapes);
    for i = 1:length(names)
        name = names{i};
        s = shapes.(name);
        if iscell(s)
            for i = 1:length(s)
                if ischar(s{i})
                    s{i} = mpctools.structget(N, s{i}, NaN());
                end
            end
            s = cell2mat(s);
        end
        out.(name) = s;
    end
end%function

