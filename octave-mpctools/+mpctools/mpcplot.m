function [vals, xax, uax] = mpcplot(varargin)
% `[vals, xax, uax] = mpcplot(x, u, [t], [xsp], [xc], [tc], ...)`
%
% Makes a plot of the given MPC solution.
%
% Arguments are as follows:
%
% * `x` : Array
%
%   > Array of size `[Nx, Nt + 1]` that gives the values of the states at the
%   discrete time points. Plotted normally.
%
% * `u` : Array
%
%   > Array of size `[Nu, Nt]` that gives the values of the controls during
%   the discrete time windows. Plotted using zero-order hold (stairstep). If
%   omitted, only x is plotted.
%
%   > Alternatively, can be size `[Nu, Nt + 1]`, which will lead to a first-
%   order hold plot.
%
% * `t` : Vector [`0:Nt`]
%
%   > Vector giving the time values at each of the discrete time points.
%
% * `xsp` : Array
%
%   > If given, defines the setpoints for the state. If size `[Nx, Nt]`, it
%   is plotted using a zero-order hold (stairstep). Otherwise, it must be
%   size `[Nx, Nt + 1]`, and it is plotted using a first-order hold (normal).
%
% * `xc` : Array
%
%   > Array of states at (interior) collocation points. Must be size
%   `[Nx, Nc, Nt]`. Used to interpolate the plot of $x$.
%
% * `tc` : Array
%
%   > Array of time points for the (interior) collocation points. Must be size
%   `[Nc, Nt]`. If not given, it is assumed that the collocation times are
%   roots of Legendre polynomials.
%
% * `usp` : Array
%
%   > Similar to `xsp`, but for `u`.
%
% * `plot` : Logical [`true`]
%
%   > Whether or not to actually perform the plotting. If `false`, this function
%   will simply return a struct with all the (properly sized) values that would
%   have been plotted.
%
% * `marker` : String [`''`]
% * `xmarker` : String
% * `umarker` : String
%
%   > Plot markers to use. The value of `marker` sets the marker to use for
%   both $x$ and $u$. To set them individually, provide `xmarker` and/or
%   `umarker`, which overrides `marker` for each variable.
%
%   > Pass the empty string (`''`) to not use any marker.
%
% * `collocmarker` : String [`''`]
%
%   > Marker to use for collocation points (if provided). Default is to use no
%   marker.
%
% * `fig` : Figure Handle
%
%   > Figure handle to use for plotting. Default is to make a new figure. This
%   is useful if you wish to plot multiple datasets on the same figure.
%
% * `title` : String [`''`]
% * `timelabel` : String [`'Time'`]
%
%   > Strings to use for the window title and bottom x-axis (time-axis) label.
%
% * `legend` : String [`''`]
%
%   > String to use as a legend entry. If nonempty, will add a legend to the
%   figure.
%
% * `legendloc` : String [`'North'`]
%
%   > String that specifies where legend should be. Should be in the format of
%   Octave/Matlab legend locations (e.g., `'West'`, `'NorthEast'`, etc.).
%
% * `color` : String or RGB array
% * `spcolor` : String or RGB array
%
%   > Colors to use for data (i.e., $x$ and $u$) and setpoint respectively.
%   Both default to black.
%
% * `xnames` : Function or Cell Array of Strings 
% * `unames` : Function or Cell Array of Strings
%
%   > Defines labels for the states and controls respectively. Used as y-axis
%   labels on each subplot.
%
%   > If a function, for each label, the function is called with a single
%   argument of the integer index of the variable. If a cell array of strings,
%   the strings are used directly.
%
% * `labelrot` : Integer
%
%   > Sets the rotation for y-axis labels (in degrees). By default, 0 rotation
%   if all labels are short (five characters or fewer), otherwise 90.
%
%
% * `linestyle` : String [`'-'`]
% * `splinestyle` : String [`':'`]
%
%   > Chooses which type of line to use for the plots. (`splinestyle` is line
%   style for setpoint `xsp` if given).
%
% The output `vals` is a struct with fields "x", "u", "t", etc. that have all
% data properly resized to have `Nt + 1` time points. It is useful if you want
% to make your own plots but don't want to have to reshape everything yourself.
%
% Other outputs `xax` and `uax` are vectors of axes handles for `x` and `u`
% respectively. Note that they will be empty if `plot` is `false`.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('x', 'required', 'numeric');
    parser.add('u', [], 'numeric');
    parser.add('t', [], {'numeric', 'row'});
    parser.add('xsp', [], 'numeric');
    parser.add('xc', [], 'numeric');
    parser.add('tc', [], 'numeric');
    parser.add('usp', [], 'numeric');
    parser.add('plot', true(), 'bool');
    parser.add('marker', '', 'str');
    parser.add('xmarker', 'yes', 'str');
    parser.add('umarker', 'yes', 'str');
    parser.add('collocmarker', '', 'str');
    parser.add('fig', []);
    parser.add('title', '', 'str');
    parser.add('timelabel', 'Time', 'str');
    parser.add('legend', '', 'str');
    parser.add('legendloc', 'North', 'str');
    parser.add('color', 'k');
    parser.add('spcolor', []);
    parser.add('xnames', [], @namechecker);
    parser.add('unames', [], @namechecker);
    parser.add('labelrot', [], {'scalar', 'numeric'});
    parser.add('linestyle', '-', 'str');
    parser.add('splinestyle', ':', 'str');
end
args = parser.parse(varargin{:});

% Get sizes.
Nx = size(args.x, 1);
Nt = size(args.x, 2) - 1;
if isempty(args.u)
    Ncols = 1;
    Nu = 0;
    uplot = @plot;
else
    Ncols = 2;
    Nu = size(args.u, 1);
    switch size(args.u, 2)
    case Nt
        uplot = @stairs;
        args.u = [args.u, NaN(size(args.u(:,end)))];
    case Nt + 1;
        uplot = @plot;
    otherwise
        error('Incorrect number of time points in u!');
    end
end
if isempty(args.t)
    args.t = 0:Nt;
end
vals = struct('x', args.x, 'u', args.u, 't', args.t); 

% Decide colors.
if isempty(args.spcolor)
    args.spcolor = args.color;
end

% Choose markers.
if isequal(args.xmarker, 'yes')
    args.xmarker = args.marker;
end
if isequal(args.umarker, 'yes')
    args.umarker = args.marker;
end

% Decide names.
if isempty(args.xnames)
    args.xnames = @(n) sprintf('x_{%d}', n);
end
if isa(args.xnames, 'function_handle')
    args.xnames = arrayfun(args.xnames, 1:Nx, 'UniformOutput', false());
elseif length(args.xnames) ~= Nx
    error('xnames must be length Nx!');
end

if isempty(args.unames)
    args.unames = @(n) sprintf('u_{%d}', n);
end
if isa(args.unames, 'function_handle')
    args.unames = arrayfun(args.unames, 1:Nu, 'UniformOutput', false());
elseif length(args.unames) ~= Nu
    error('unames must be length Nu!');
end

% Decide label rotation.
if isempty(args.labelrot)
    if max(labellength([args.xnames, args.unames])) <= 5
        args.labelrot = 0;
    else
        args.labelrot = 90;
    end
end

% Check collocation.
if ~isempty(args.xc)
    usecollocation = true();
    [X, T] = mpctools.smushcolloc(args.x, args.xc, args.t, args.tc);
    Nc = size(args.xc, 2);
    vals.X = X;
    vals.T = T;
    vals.xc = X;
    vals.xc(:,1:(Nc + 1):end) = [];
    vals.tc = T;
    vals.tc(1:(Nc + 1):end) = [];
else
    usecollocation = false();
end

% Check setpoint.
if ~isempty(args.xsp)
    if isequal(size(args.xsp), [Nx, Nt])
        vals.xsp = [args.xsp, args.xsp(:,end)];
        vals.plotxsp = @stairs;
    elseif isequal(size(args.xsp), [Nx, Nt + 1])
        vals.xsp = args.xsp;
        vals.plotxsp = @plot;
    else
        error('Invalid size for xsp!');
    end
else
    vals.xsp = [];
end
if ~isempty(args.usp)
    if isequal(size(args.usp), [Nu, Nt])
        vals.usp = [args.usp, args.usp(:,end)];
        vals.plotusp = @stairs;
    elseif isequal(size(args.usp), [Nu, Nt + 1])
        vals.usp = args.usp;
        vals.plotusp = @plot;
    else
        error('Invalid size for usp!');
    end
else
    vals.usp = [];
end

% Make plot.
xax = cell(Nx, 1);
uax = cell(Nu, 1);
if args.plot
    if isempty(args.fig)
        args.fig = figure();
    end
    set(0, 'CurrentFigure', args.fig); % figure() overrides visibility.
    vals.fig = args.fig;
    Nrows = max(Nx, Nu);
    
    % Do x.
    for i = 1:Nx
        xax{i} = subplotij(Nrows, Ncols, i, 1);
        if usecollocation
            plot(vals.T, vals.X(i,:), args.linestyle, 'color', args.color);
            plot(vals.t, vals.x(i,:), args.xmarker, 'color', args.color);
            if ~isempty(args.collocmarker)
                plot(vals.tc, vals.xc(i,:), args.collocmarker, ...
                     'color', args.color);
            end
        else
            plot(vals.t, vals.x(i,:), [args.linestyle, args.xmarker], ...
                 'color', args.color);
        end
        if ~isempty(vals.xsp)
            vals.plotxsp(vals.t, vals.xsp(i,:), args.splinestyle, ...
                         'color', args.spcolor);
        end
        ylabel(args.xnames{i}, 'rotation', args.labelrot);
    end
    if ~isempty(args.timelabel)
        xlabel(args.timelabel);
    end
    
    % Do u.
    for i = 1:Nu
        uax{i} = subplotij(Nrows, Ncols, i, 2);
        uplot(vals.t, vals.u(i,:), args.linestyle, 'color', args.color);
        if ~isempty(args.umarker)
            plot(vals.t, vals.u(i,:), args.umarker, 'color', args.color);
        end
        if ~isempty(vals.usp)
            vals.plotusp(vals.t, vals.usp(i,:), args.splinestyle, ...
                         'color', args.spcolor);
        end
        ylabel(args.unames{i}, 'rotation', args.labelrot);
    end
    if ~isempty(args.timelabel)
        xlabel(args.timelabel);
    end
    
    % Add title.
    if ~isempty(args.title)
        mpctools.windowtitle(args.title);
    end
    
    % Add legend.
    if ~isempty(args.legend)
        if Ncols == 1
            i = Nx;
            j = 1;
        elseif Nu > Nx
            i = Nx + 1;
            j = 1;
        elseif Nx > Nu
            i = Nu + 1;
            j = 2;
        else
            i = Nu;
            j = 2;
        end
        subplotij(Nrows, Ncols, i, j);
        plot(NaN(), NaN(), [args.linestyle, args.xmarker], ...
             'color', args.color, 'displayname', args.legend);
        legend('show');
        legend('location', args.legendloc);
        if Ncols == 2 && Nx ~= Nu
            axis('off');
        end
    end
end
xax = vertcat(xax{:});
uax = vertcat(uax{:});

end%function

function ax = subplotij(Nrows, Ncols, i, j)
    % A sensibly-defined subplot function.
    p = (i - 1)*Ncols + j;
    ax = subplot(Nrows, Ncols, p);
    hold('on');
end%function

function names = namechecker(names)
    % Checks that names is a cell array of strings or a single function handle.
    if isa(names, 'function_handle')
        % pass
    elseif iscell(names)
        isstrs = cellfun(@ischar, names);
        if ~all(isstrs(:))
            error('All cell entries must be strings!');
        end
        names = names(:)';
    else
        error('Must be a cell array of strings or a function handle!');
    end
end%function

function lengths = labellength(labels)
    % Returns an array the same size as cell array labels with estimates of
    % how long the printed strings will be.
    %
    % The function attempts to account for interpreted sequences whose raw
    % character counts are larger than how many characters will actually be
    % printed (e.g., subscripts with braces), but it is not exact.
    lengths = zeros(numel(labels), 1);
    for i = 1:numel(labels)
        label = regexprep(labels{i}, '[_^]\{([^}]*)\}', '$1');
        label = regexprep(label, '\\[A-Za-z]+\s*', '#');
        lengths(i) = length(label);
    end
end%function

