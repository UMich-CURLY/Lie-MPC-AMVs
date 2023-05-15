function [sym, varlist, parlist] = getocpsymbols(casaditype, shapes, N, ...
                                                 defaultvar, customvar, ...
                                                 parval, finalx, finaly)
% [sym, varnames, parnames] = getocpsymbols(shapes, N, defaultvar, customvar,
%                                           parval, [finalx=true], [finaly=true])
%
% Returns symbolic struct of variables and parameters.
%
% Remaining arguments are cell arrays giving the names of variables and
% parameters.
narginchk(6, 8);
if ~all(cellfun(@isstruct, {shapes, N, parval}))
    error('shapes, N, and parval must be structs!');
elseif ~iscell(defaultvar) || ~iscell(customvar)
    error('defaultvar and customvar must be cell arrays!');
end
if nargin() < 7
    finalx = true();
end
if nargin() < 8
    finaly = true();
end

% Sort out default variables. Note that we have to check consistency with the
% shapes that were inferred from the functions.
[varshapes, varrepeats, defaultshapes] = ...
    mpctools.generalvariableshapes_(defaultvar, N, finalx, finaly);
cellfun(@(v) checkconsistency(v, shapes, varshapes), defaultvar);

% Sort out custom variables.
badcustom = intersect(defaultvar, customvar);
if ~isempty(badcustom)
    error('Custom variables {%s} clash with default variables!', ...
          mpctools.row2str(badcustom));
end
for i = 1:length(customvar)
    v = customvar{i};
    if ~isfield(shapes, v)
        error('Custom variable "%s" never used! Unable to determine size!', v);
    end
    varshapes.(v) = shapes.(v);
    varrepeats.(v) = zeros(0, 1);
end
varlist = [defaultvar(:)', customvar(:)'];

% Sort out parameters. Shapes may be in the varshapes struct.
if ismember('Du', varlist)
    if ~isfield(parval, 'uprev')
        parval.uprev = zeros(varshapes.Du);
    end
    if ~isfield(shapes, 'uprev')
        shapes.uprev = varshapes.Du;
    end
end
parlist = union(fieldnames(parval), setdiff(fieldnames(shapes), varlist));
keep = true(size(parlist));
parshapes = struct();
for i = 1:length(parlist)
    p = parlist{i};
    if isfield(shapes, p)
        parshapes.(p) = shapes.(p);
    elseif isfield(defaultshapes, p) && ~any(isnan(defaultshapes.(p)))
        parshapes.(p) = defaultshapes.(p);
        checkconsistency(p, parshapes, shapes);
    else
        keep(i) = false();
        warning('Skipping unused parameter "%s".', p);
    end
end
parlist = parlist(keep);
parrepeats = calcparrepeats_(parshapes, parval);

% Check for overlaps.
overlaps = intersect(varlist, parlist);
if ~isempty(overlaps)
    error('The following parameter names clash with variable names: {%s}!', ...
          mpctools.row2str(overlaps));
end

% Get combined struct.
shapes = mpctools.structupdate(varshapes, parshapes);
repeats = mpctools.structupdate(varrepeats, parrepeats);
sym = mpctools.getSymStruct(casaditype, shapes, repeats);

end%function

function repeats = calcparrepeats_(shapes, vals)
    % repeats = calcparrepeats_(shapes, vals)
    %
    % Calculates the appropriate number of repeats for a struct of parameters.
    %
    % shapes should be a struct of parameter names and sizes for the individual
    % symbolic units (e.g., size of vector or matrix). vals should be a struct
    % of values as an ND array (with symbols concatenated first, and then
    % repeats).
    narginchk(2, 2);
    if ~isstruct(shapes) || ~isstruct(vals)
        error('shapes and vals must both be structs!');
    end

    names = fieldnames(shapes);
    given = true(size(names));
    repeats = struct();
    
    for i = 1:length(names)
        n = names{i};
        vecsize = shapes.(n);
        given(i) = isfield(vals, n);
        if given(i)
            v = vals.(n);
            [repeats.(n), errmsg] = calcrepeat(vecsize, size(v));
            if ~isempty(errmsg)
                error('Incorrect size for parameter "%s": %s\n', n, errmsg);
            end
        end
    end
    
    % Print error message if any parameters not given.
    missing = names(~given);
    if ~isempty(missing)
        error('Some parameters were used but no value was given: {%s}', ...
              mpctools.row2str(missing));
    end
end%function

function okay = checkconsistency(v, struct1, struct2)
    % Checks that struct1.(v) and struct2.(v) are equal. If either are empty,
    % doesn't do anything.
    if isfield(struct1, v) && isfield(struct2, v);
        check1 = struct1.(v);
        check2 = struct2.(v);
        okay = isequal(check1, check2);
        if ~okay
            error('Size mismatch for "%s": expected [%s], got [%s]!', v, ...
                  mpctools.row2str(check1), mpctools.row2str(check2));
        end
    else
        okay = true();
    end
end%function

function [repeat, errmsg] = calcrepeat(vecsize, fullsize)
    % Calculates the offset index in fullsize that defines the repeats.
    errmsg = '';
    if isequal(vecsize, [1, 1])
        % Special case because this is only an error if given something empty.
        if isequal(fullsize, [1, 1])
            repeat = zeros(1, 0);
        elseif length(fullsize) == 2 && fullsize(2) == 1
            repeat = fullsize(1);
        elseif prod(fullsize) > 0
            repeat = fullsize;
        else
            errmsg = sprintf('given empty array with size [%s]!', ...
                             mpctools.row2str(fullsize));
            repeat = [];
        end
    else
        % Here we have to guess whether singleton dimensions are significant.
        expectsize = [vecsize(1), NaN()];
        if vecsize(2) ~= 1 || fullsize(2) == 1
            expectsize(2) = vecsize(2);
            offset = 2;
        else
            offset = 1;
        end
        s1 = expectsize(1:offset);
        s2 = fullsize(1:offset);
        if isequal(s1, s2)
            repeat = fullsize((offset + 1):end);
        else
            errmsg = sprintf('expected size [%s, ...]; got [%s, ...].', ...
                             mpctools.row2str(s1), mpctools.row2str(s2));
            repeat = [];
        end
    end
end%function

