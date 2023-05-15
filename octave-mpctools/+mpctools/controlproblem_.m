function solver = controlproblem_(varargin)
% solver = controlproblem_(...)
%
% Builds a general optimal control problem.
%
% This function is private and should not be invoked directly. Use `nmpc()`,
% `nmhe()`, `sstarg()`, or `parest()` instead.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('N', 'required', 'struct');
    parser.add('var', 'required', 'cell');
    parser.add('sym', 'required', 'struct');
    parser.add('funcs', 'required', 'struct');
    parser.add('funcargs', 'required', 'struct');
    parser.add('par', {}, 'cell');
    parser.add('parval', [], 'struct');
    parser.add('lb', struct(), 'struct');
    parser.add('ub', struct(), 'struct');
    parser.add('guess', struct(), 'struct');
    parser.add('obj', 0);
    parser.add('con', {}, 'cell');
    parser.add('Delta', NaN(), {'scalar', 'pos'});
    parser.add('discreteness', struct(), 'struct');
    parser.add('discretel', true(), 'bool');
    parser.add('infercolloc', {'guess', 'lb', 'ub'}, 'cell');
    parser.add('singleshooting', false(), 'bool');
    parser.add('solvertype', 'default', 'str');
    parser.add('finalstagecost', false(), 'bool');
    parser.add('collocjump', struct(), 'struct');
    parser.add('**kwargs');
end
[args, extraargs] = parser.parse(varargin{:});

% Get symbolic variables.
sym = args.sym;
var = mpctools.structupdate(struct(), sym, args.var);
par = mpctools.structupdate(struct(), sym, args.par);
parval = args.parval;

% Decide whether we're using collocation.
N = args.N;
usecollocation = isfield(N, 'c') && (N.c > 0);

weights = struct();
if usecollocation
    if isnan(args.Delta)
        error('Argument "Delta" must be given to use collocation!');
    end
    [weights.r, weights.A, ~, weights.q] ...
        = mpctools.collocweights(N.c,  'right', 'left');
        
    % Infer values for collocation points.
    for i = 1:length(args.infercolloc)
        field = args.infercolloc{i};
        if ~ismember(field, {'guess', 'lb', 'ub'})
            error('Invalid value for infercolloc');
        end
        for v = ['x', 'z']
            if isfield(args.(field), v) && ~isfield(args.(field), [v, 'c'])
                args.(field).([v, 'c']) = mpctools.infercolloc(weights.r, ...
                                                               args.(field).(v));
            end
        end
    end
else
    weights.r = [];
end
if isnan(args.Delta)
    args.Delta = 1;
end

% Get default values for bounds and update from user arguments.
varndims = struct('xc', 2, 'zc', 2);

lb = mpctools.varstruct2numstruct(var, -inf(), varndims);
lb = updatebounds(lb, args.lb);

ub = mpctools.varstruct2numstruct(var, inf(), varndims);
ub = updatebounds(ub, args.ub);

guess = mpctools.varstruct2numstruct(var, 0, varndims);
guess = updatebounds(guess, args.guess);

discreteness = mpctools.varstruct2numstruct(var, false(), varndims);
discreteness = updatebounds(discreteness, args.discreteness);

% Get function arguments and expand symbols for constraints.
funcargs = args.funcargs;
if usecollocation
    symc = mpctools.collocmerge_(sym, funcargs, args.collocjump);
end

% Build model constraints.
if isfield(args.funcs, 'f') && ~isempty(args.funcs.f)
    f = args.funcs.f;
    if usecollocation
        if args.singleshooting
            error('singleshooting is not available with collocation!');
        end
        model = mpctools.collocconstraints(f, symc, funcargs.f, ...
                                           weights.A, args.Delta);
    else
        if args.singleshooting
            model = struct('con', {var.x(:,2:end)}, ...
                           'lb', lb.x(:,2:end), ...
                           'ub', ub.x(:,2:end));
            var.x = var.x(1);
            lb.x = lb.x(:,1);
            ub.x = ub.x(:,1);
            guess.x = guess.x(:,1);
            discreteness.x = discreteness.x(:,1);
        else
            model = mpctools.multipleshooting(f, sym, funcargs.f);
        end
    end
else
    model = mpctools.emptycon();
end

% Algebraic state constraints.
if isfield(args.funcs, 'g') && ~isempty(args.funcs.g)
    g = args.funcs.g;
    if usecollocation
        [algebra, collocalg] = mpctools.algebraconstraints(g, sym, ...
                                                           funcargs.g, symc);
    else
        algebra = mpctools.algebraconstraints(g, sym, funcargs.g);
        collocalg = mpctools.emptycon();
    end
else
    algebra = mpctools.emptycon();
    collocalg = mpctools.emptycon();
end

% Measurement constraints.
if isfield(args.funcs, 'h') && ~isempty(args.funcs.h)
    measurements = mpctools.measurementconstraints(args.funcs.h, sym, ...
                                                   funcargs.h);
else
    measurements = mpctools.emptycon();
end

% Path constraints.
if isfield(args.funcs, 'e') && ~isempty(args.funcs.e)
    pathcon = mpctools.pathconstraints(N.t, args.funcs.e, sym, funcargs.e);
else
    pathcon = mpctools.emptycon();
end

% Rate-of-change constraints.
deltacon = {};
for v = ['x', 'u', 'd']
    if isfield(var, ['D', v])
        deltacon = [deltacon, {['delta', v], mpctools.rocconstraints(sym, v)}];
    end
end

% Absolute value constraints and bounds.
needabs = intersect({'absx', 'absxc', 'absz', 'abszc', 'absu', 'absDx', ...
                     'absDz', 'absDu', 'absw', 'absv', 'absy', ...
                     'absd', 'absDd'}, fieldnames(var));
abscon = cell(1, 2*length(needabs));
for i = 1:length(needabs)
    v = needabs{i};
    abscon{2*i - 1} = v;
    abscon{2*i} = mpctools.absconstraints(var, v(4:end));
    lb.(v) = zeros(size(lb.(v)));
    ub.(v) = inf(size(ub.(v)));
end

% Assemble all constraints and constraint bounds.
allcon = struct('model', model, deltacon{:}, ...
                'path', pathcon, 'measurements', measurements, ...
                'algebra', algebra, 'collocalgebra', collocalg, ...
                abscon{:}, args.con{:});
[con, conlb, conub] = mpctools.constraintmerge(allcon);

% Build objective function.
if isfield(args.funcs, 'l') && ~isempty(args.funcs.l)
    l = args.funcs.l;
    if usecollocation && ~args.discretel
        obj = mpctools.quadstagecosts(N.t, l, symc, funcargs.l, ...
                                      args.obj, weights.q, args.Delta);
    elseif args.discretel
        obj = mpctools.sumstagecosts(N.t, l, sym, funcargs.l, args.obj);
    else
        error('discretel=False requires N.c > 0!');
    end
else
    obj = 0;
end

% Create ControlSolver object.
kwargs = struct('horizon', N.t, 'var', var, 'varlb', lb, 'varub', ub, ...
                'varguess', guess, 'obj', obj, 'con', con, ...
                'conlb', conlb, 'conub', conub, 'collocpts', weights.r, ...
                'Delta', args.Delta, 'parval', parval, 'par', par, ...
                'discreteness', discreteness, ...
                'singleshooting', args.singleshooting);
switch args.solvertype
case 'mhe'
    solver = mpctools.MHESolver('**', kwargs, '**', extraargs);
otherwise
    solver = mpctools.ControlSolver('**', kwargs, '**', extraargs);
end

end%function

function s = updatebounds(s, bounds)
    % s = updatebounds(s, bounds)
    %
    % Updates the bounds in s using the fields of bounds.
    fields = intersect(fieldnames(s), fieldnames(bounds));
    for i = 1:length(fields)
        f = fields{i};
        datasize = size(bounds.(f));
        defaultsize = size(s.(f));
        [okay, rep] = mpctools.isbsxable(datasize, defaultsize);
        if ~okay
            error('Unknown size for field "%s" of %s!', f, inputname(1));
        end
        s.(f) = repmat(bounds.(f), rep);
    end
end%function

