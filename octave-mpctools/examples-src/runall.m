function failures = runall(doplots)
% `failures = runall([showplots=False])`
%
% Runs all of the example scripts in this directory and returns a struct of
% errors for any scripts that failed.
%
% If `showplots` is True, plots will be displayed as they are created (note that
% they will steal focus, which can be annoying). Otherwise, they will remain
% hidden until the script finishes (at which point they will also steal focus).
% Finally, if `showplots` is the string 'off', then they will never be shown.
narginchk(0, 1);
clearplots = false();
if nargin() < 1
    doplots = false();
elseif isequal(doplots, 'off')
    doplots = false();
    clearplots = true();
end

% Set some local printing options if octave.
if mpctools.isOctave()
    more('off');
    page_output_immediately(true(), 'local');
    page_screen_output(false(), 'local');
end

% Start with a dummy NLP to display the IPOPT splash message.
x = casadi.SX.sym('x');
nlp = struct('x', x, 'f', x.^2);
ipopt = casadi.nlpsol('ipopt', 'ipopt', nlp, ...
                      struct('print_time', false(), ...
                             'ipopt', struct('print_level', 0)));
ipopt('x0', 0);

% Switch off plotting.
DefaultFigureVisible = get(0, 'DefaultFigureVisible');
if ~doplots
    set(0, 'DefaultFigureVisible', 'off'); % Don't show figures.
end
initialfigures = get(0, 'children');

% Now actually run the examples.
examples = {
    'ballmaze';
    'blockstructure';
    'cstr';
    'cstr_startup';
    'cstr_nmpc_nmhe';
    'demcharge';
    'econmpc';
    % 'fishing'; % Exclude for now until Bonmin output can be suppressed.
    'fullinformationexample';
    'mheexample';
    'nmheexample';
    'nmpcexample';
    'onenorm';
    'rocpenalty';
    'simplemhe';
    'simplempc';
    'softconstraints';
    % 'sstargexample'; % Exclude for now until Bonmin output can be suppressed.
    'timevaryingmpc';
    'vdposcillator';
    'rounding';
    'gainscheduling';
    'daeexample';
    'springpendulummhe';
    'timevaryingmhe';
    'priorupdates';
    'predatorprey';
    'mhe_discretization';
};
fprintf('\nRunning MPCTools example scripts.\n\n');
failures = struct();
if mpctools.isOctave()
    runscript = @runscript_octave;
else
    runscript = @runscript_matlab;
end
for i = 1:length(examples)
    ex = examples{i};
    fprintf('* %s ... ', ex);
    extime = tic();
    err = runscript(ex, doplots);
    extime = toc(extime);
    if isempty(err)
        fprintf('success');
    else
        fprintf('FAILURE');
        failures.(ex) = err;
    end
    fprintf(' (%.4g s)\n', extime);
end

% Make figures show up.
if ~doplots
    newplots = setdiff(get(0, 'children'), initialfigures);
    for i = 1:length(newplots)
        if clearplots
            close(newplots(i));
        else
            set(newplots(i), 'Visible', 'on');
        end
    end
end

% Switch plots back on.
set(0, 'DefaultFigureVisible', DefaultFigureVisible);

end%function

function err = runscript_matlab(script, ~)
    % err = runscript_matlab(script, ~)
    %
    % Attempts to run a script without displaying anything and returns any error
    % message.
    try
        evalc(script);
        err = [];
    catch err
        % Pass.
    end
end%function

function err = runscript_octave(script, doplots)
    % err = runscript_matlab(script, doplots)
    %
    % Attempts to run a script and returns any error message.
    %
    % Most output should be suppressed, although we can't hit everything.
    noop = str2func('mpctools.noop');
    fprintf = noop;
    disp = noop;
    display = noop;
    if ~doplots
        % Monkey-patch the figure command so plots don't show up.
        figure = @(varargin) builtin('figure', varargin{:}, 'visible', 'off');
    end
    verb = mpctools.MAX_VERBOSITY(0); % Disable solver output.
    try
        eval(script);
        err = [];
    catch err
        % Pass.
    end
    mpctools.MAX_VERBOSITY(verb); % Reset verbosity.
end%function

