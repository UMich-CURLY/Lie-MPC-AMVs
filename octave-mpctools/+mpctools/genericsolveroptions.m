function opts = genericsolveroptions(varargin)
% opts = genericsolveroptions(solver, timelimit, verbosity, isQP, maxiter)
%
% Returns a struct of options based on the given solver and generic names of
% this function. This allows solver-invariant names to be used for common
% options.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('solver');
    parser.add('timelimit', [], {'scalar', 'pos'});
    parser.add('verbosity', [], {'int', 'nonneg'});
    parser.add('isQP', [], 'bool');
    parser.add('maxiter', [], {'int', 'pos'});
    parser.add('discretevar', [], {'logical', 'col'});
end
args = parser.parse(varargin{:});

% Set option names.
sopts = struct();
opts = struct();
switch args.solver
case {'ipopt', 'bonmin'}
    sopts.max_cpu_time = args.timelimit;
    sopts.print_level = args.verbosity;
    if isequal(args.solver, 'bonmin')
        if ~isempty(args.verbosity)
            % Note that these don't seem to have any effect, probably due to an
            % issue with Casadi's bonmin interface.
            sopts.bb_log_level = floor(args.verbosity*5/12);
            sopts.fp_log_level = floor(args.verbosity*2/12);
            sopts.oa_log_level = floor(args.verbosity*2/12);
            sopts.nlp_log_level = floor(args.verbosity*2/12);
        end
        sopts.time_limit = args.timelimit;
        sopts.iteration_limit = args.maxiter;
    else
        sopts.max_iter = args.maxiter;
    end
    if ~isempty(args.isQP) && args.isQP
        qpoptions = {'jac_c_constant', 'jac_d_constant', 'hessian_constant'};
        for i = 1:length(qpoptions)
            sopts.(qpoptions{i}) = 'yes';
        end
    end
    if ~isempty(args.verbosity)
        opts.print_time = (args.verbosity > 2);
    end
    opts.eval_errors_fatal = true();
case 'qpoases'
    opts.CPUtime = args.timelimit;
    if ~isempty(args.verbosity)
        if args.verbosity >= 9
            plevel = 'high';
        elseif args.verbosity >= 6
            plevel = 'medium';
        elseif args.verbosity >= 3
            plevel = 'low';
        else
            plevel = 'none';
        end
        opts.printLevel = plevel;
    end
    if ~isempty(args.maxiter)
        warning('qpoases does not support maxiter!');
    end
    opts.sparse = true();
    sopts = []; % qpoases doesn't use this one.
case 'blocksqp'
    % TODO: fix blocksqp options.
    opts.time = args.timelimit;
    opts.max_iter = args.maxiter;
    sopts = [];
case 'gurobi'
    % TODO: fix this when an options interface is added to Casadi.
    warning('Gurobi options are not yet available!');
    %opts.TimeLimit = args.timelimit;
    %opts.OutputFlag = ceil(args.verbosity/12);
    sopts = []; % gurobi doesn't use these.
otherwise
    warning('Solver %s is not currently supported!', args.solver);
end

% Get rid of any that are the empty matrix.
opts = purgeempty(opts);
if ~isempty(sopts)
    sopts = purgeempty(sopts);
    opts.(args.solver) = sopts;
end

% Check for discrete variables.
if ~isempty(args.discretevar)
    if ismember(args.solver, {'bonmin', 'gurobi', 'cplex'})
        opts.discrete = args.discretevar;
    elseif any(args.discretevar)
        error('Solver %s does not support discrete variables!', args.solver);
    end
end

end%function

function s = purgeempty(s)
    % Removes the empty fields of s.
    fields = fieldnames(s);
    remove = false(size(fields));
    for i = 1:length(fields)
        f = fields{i};
        if isempty(s.(f))
            remove(i) = true();
        end
    end
    s = rmfield(s, fields(remove));
end%function

