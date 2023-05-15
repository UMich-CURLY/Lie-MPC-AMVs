function varargout = getLinearizedModel(varargin)
% `model = getLinearizedModel(f, args, names, [Delta], [deal=false])`
%
% Linearizes the model
%
% > $\frac{dx}{dt} = f(x, u, ...)$
%
% the point `{xss, uss, ...}` given in `args`. `names` should be a cell array to
% define the fields into which each matrix should go. For example, to linearize
% the model
%
% > $f(x,u) \approx Ax + Bu$
%
% you should use `names = {'A', 'B'}`, and `model` will be a struct with
% fields "A" and "B".
%
% `f` should be either a casadi.Function object or a native function handle.
%
% If `Delta` is given, the model is also discretized with that timestep.
%
% By default, the return value model is a struct whose fields are the strings
% in `names`. However, if `deal` is set to `true`, the function will return
% the individual matrices as multiple outputs, e.g.
% ```
%     [A, B] = getLineraizedModel(f, {xss, uss}, {'A', 'B'}, 'deal', true())
% ```
% Sometimes, this syntax is more convenient.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('f', 'required', 'func');
    parser.add('args', 'required', {'cell', 'row'});
    parser.add('names', [], 'cell');
    parser.add('Delta', [], {'scalar', 'pos'});
    parser.add('deal', false(), 'bool');
end
args = parser.parse(varargin{:});

% Get names.
nargs = length(args.args);
if nargs == 0
    error('f must take at least one argument!');
end
if isempty(args.names)
    if nargs == 2
        args.names = {'A', 'B'};
    else
        args.names = cell(1, nargs);
        args.names{1} = 'A';
        for i = 1:(nargs - 1)
            args.names{i + 1} = sprintf('B_%d', i);
        end
    end
end

% Check if f is a casadi function. If not, need to make it one.
if ~mpctools.iscasadifunc(args.f)
    sizes = cellfun(@numel, args.args);
    args.f = mpctools.getCasadiFunc(args.f, sizes);
end

% Get jacobians.
jblocks = cellfun(@(n) sprintf('jac:%s:%s', args.f.name_out(0), n), ...
                  mpctools.names2cell_(args.f.name_in()), ...
                  'UniformOutput', false());
jac = args.f.factory('jac_f', args.f.name_in, jblocks);
jacobians = cell(1, nargs);
[jacobians{:}] = jac(args.args{:});
for i = 1:nargs
    jacobians{i} = full(jacobians{i});
end

% Discretize if requested.
if ~isempty(args.Delta)
    mats = mpctools.c2d(args.Delta, jacobians{1}, eye(size(jacobians{1})));
    jacobians{1} = mats.A;
    for i = 2:nargs
        jacobians{i} = mats.B*jacobians{i};
    end
end

% Choose return value.
if args.deal
    varargout = jacobians;
else
    model = struct();
    for i = 1:nargs
        model.(args.names{i}) = jacobians{i};
    end
    varargout = {model};
end

end%function
