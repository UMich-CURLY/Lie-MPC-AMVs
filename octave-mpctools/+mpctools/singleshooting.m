function vars = singleshooting(vars, f, fargs, h, hargs, trange)
% vars = singleshooting(vars, f, fargs, [h], [hargs], [trange])
%
% Replaces the "x" field in vars with the corresponding single-shooting
% expressions. If h and hargs are also given, the "y" field in vars is also
% updated appropriately.
%
% Note that vars.x needs to be the correct length, e.g., by first initializing
% it for multiple shooting.
narginchk(3, 6);
needh = nargin() > 3;
if needh
    if nargin() < 5 || isempty(hargs)
        error('If h is given, hargs must also be given!');
    end
end
if ~isfield(vars, 'x')
    error('Field "x" is missing from vars!');
end
if nargin() < 6
    trange = 1:length(vars.x);
end

% Compute expressions for x variables.
for i = trange(1:(end - 1))
    args = mpctools.getargs_(i, vars, fargs);
    try
        vars.x{i + 1} = f(args{:});
    catch err
        error('evaluating f(%s): %s', mpctools.row2str(fargs), err.message);
    end
end

% Compute expressions for y variables.
if needh
    if ~isfield(vars, 'y')
        error('Field "y" is missing from vars!');
    end
    for i = trange
        args = mpctools.getargs_(i, vars, hargs);
        try
            vars.y{i} = h(args{:});
        catch err
            error('evaluating h(%s): %s', mpctools.row2str(hargs), ...
                  err.message);
        end
    end
end

end%function

