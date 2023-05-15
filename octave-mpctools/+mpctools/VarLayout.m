classdef VarLayout < handle
% Class for holding a layout of variables. Contains information about how to map
% from the (possibly multidimensional) variable space into the single large
% vector that is needed by the solvers.

properties (GetAccess = public, SetAccess = private)
    names; % Cell array of names.
    
    sizes; % Cell array of sizes.
    squeezedsizes; % Cell array of sizes with singletons removed.
    Nts; % Vector of dimensions.
    maxNt; % Scalar.
    ismat; % Vector of bool.
    isconst; % Vector of bool.
    
    allinds; % Vector of integers.
    tinds; % Vector of integers.
    constinds; % Vector of integers.
    
    Nnames; % Constant. Number of different variable names.
    Nvar; % Constant. Total number of scalar variable entries.
    
    map; % Cell of indices.
    
    % Remaining parameters set by constructor.
    always_checksize;
    parameters_at_end;
    allow_squeezed;
end%properties

methods (Access = public)
    function self = VarLayout(varargin)
        % layout = VarLayout(syms, ...)
        %
        % Initializes the layout of variables.
        %
        % In the first signature, syms should be a struct with each entry either
        % a single Casadi variable (for time-invariant variables) or a cell
        % array of Casadi variables (for time-varying variables).
        %
        % Additional arguments (should be passed as 'key', value pairs) are as
        % follows:
        %
        % - 'always_checksize': when vectorizing, whether to always check the
        %                       size of each field in the struct; if false, only
        %                       the total number of elements is checked
        %                       (default False)
        %
        % - 'parameters_at_end': whether to put parameters (i.e., variables that
        %                        are time-invariant) at the end of the vector
        %                        (default True)
        %
        % - 'allow_squeezed': whether to allow numeric structs to have
        %                     singleton dimensions removed (default True)
        %
        % - 'matrixvar': cell array of variables that are matrices; used to
        %                force column vectors to be treated as matrices
        %                (default {})
        %
        % - 'ignore': cell array of fields to ignore (default {})
        %
        % - 'squeeze_row_cells': whether to squeeze the first dimension out of
        %                        row cell arrays (default True)
        %
        % Useful functions are self.vec (self.symvec) to convert from a
        % numeric (symbolic) struct to a single vector, and self.ivec to convert
        % from a vector back to a symbolic struct.
        persistent parser
        if isempty(parser)
            parser = mpctools.ArgumentParser();
            parser.add('syms', 'required', {'scalar', 'struct'});
            parser.add('parameters_at_end', true(), 'bool');
            parser.add('allow_squeezed', true(), 'bool');
            parser.add('matrixvar', {}, {'cell', 'row'});
            parser.add('ignore', {}, {'cell', 'row'});
            parser.add('squeeze_row_cells', true(), 'bool');
        end
        args = parser.parse(varargin{:});
        
        % Grab arguments.
        syms = args.syms;
        self.parameters_at_end = args.parameters_at_end;
        self.allow_squeezed = args.allow_squeezed;
        
        % Get list of names. Note that some may be ignored.
        names_ = self.sortnames(fieldnames(syms));
        keep = true(size(names_));
        for i = 1:length(names_)
            keep(i) = ~ismember(names_{i}, args.ignore);
        end
        self.names = names_(keep);
        self.Nnames = length(self.names);
        self.allinds = 1:self.Nnames;
        
        % Decide whether matrices or vectors.
        self.isconst = false(self.Nnames, 1);
        self.ismat = false(self.Nnames, 1);
        for i = self.allinds
            name = self.names{i};
            thissym = syms.(name);
            self.isconst(i) = ~iscell(thissym);
            if ~self.isconst(i)
                if isempty(thissym)
                    error('Variable "%s" is empty!', name);
                end
                thissym = thissym{1}; % Grab first element.
            end
            if ~isa(thissym, 'casadi.SX') && ~isa(thissym, 'casadi.MX')
                error('Field "%s" does not contain symbolic variables!', name);
            end
            self.ismat(i) = size(thissym, 2) > 1 ...
                            || ismember(name, args.matrixvar);
        end
        self.constinds = self.allinds(self.isconst);
        self.tinds = self.allinds(~self.isconst);
        
        % Decide sizes.
        self.sizes = cell(self.Nnames, 1);
        self.squeezedsizes = cell(self.Nnames, 1);
        for i = self.allinds
            name = self.names{i};
            thissym = syms.(name);
            if self.isconst(i)
                self.sizes{i} = size(thissym);
            else
                if self.ismat(i)
                    varsize = size(thissym{1});
                else
                    varsize = size(thissym{1}, 1);
                end
                cellsize = size(thissym);
                if args.squeeze_row_cells && isrow(thissym)
                    cellsize = cellsize(2);
                end
                totalsize = [varsize, cellsize];
                while length(totalsize) > 2 && totalsize(end) == 1
                    totalsize = totalsize(1:(end - 1));
                end
                self.sizes{i} = totalsize;
                self.squeezedsizes{i} = totalsize(totalsize ~= 1);
            end
        end
        
        % Pass things to calcmap.
        self.calcmap();
    end%function
    
    function v = vec(self, s)
        % v = self.vec(s)
        %
        % Collapses the struct s (which must be of the right shape) into a
        % single large vector. Each element of s must be a single large array.
        narginchk(2, 2);
        if self.always_checksize
            numtype = self.checknumeric(s);
        else
            numtype = class(s.(self.names{1}));
        end
        try
            v = self.vec_numeric(s, numtype);
        catch err
            % Check sizes if we didn't already do that.
            if ~self.always_checksize
                self.checknumeric(s);
            end
            if mpctools.isOctave()
                err.message = ['Error in vec: ', err.message];
            end
            rethrow(err);
        end
    end%function
    
    function v = symvec(self, s)
        % v = self.symvec(s)
        %
        % Collapses the struct s (which must be of the right shape) containing
        % Casadi symbolic variables (or cell arrays of symbolic variables) into
        % a single Casadi symbolic vector.
        narginchk(2, 2);
        self.checksize(s); % Make sure variables are the right size.
        v = self.vec_symbolic(s);
    end%function
    
    function s = ivec(self, v)
        % s = self.ivec(v)
        %
        % Returns a struct representation of the long vector v. v must be
        % numeric.
        narginchk(2, 2);
        if ~self.isarray(v)
            error('Only numeric vectors are allowed!');
        elseif numel(v) ~= self.Nvar
            error('Expected %d elements in v but got %d!', self.Nvar, numel(v));
        end
        s = struct();
        for i = self.allinds
            name = self.names{i};
            thismap = self.map{i};
            s.(name) = reshape(v(thismap), size(thismap));
        end
    end%function
    
    function m = getmap(self)
        % m = self.getmap()
        %
        % Returns the mapping for all variables.
        m = struct()
        for i = 1:self.allinds
            name = self.names{i};
            m.(name) = self.map{i};
        end
    end%function
end%methods

methods (Access = private)
    function calcmap(self)
        % Recalculates the map for the current sizes. Time-varying variables
        % are at the top, followed by time-invariant variables.
        self.Nvar = sum(cellfun(@prod, self.sizes));
        v = 1:self.Nvar;
        
        % Decide time dimensions for time-varying variables.
        self.Nts = zeros(self.Nnames, 1);
        tdims = zeros(self.Nnames, 1);
        for i = self.tinds
            s = self.sizes{i};
            tdims(i) = max(length(s), 2 + double(self.ismat(i)));
            if tdims(i) > length(s)
                self.Nts(i) = 1; % Trailing singleton dimensions.
            else
                self.Nts(i) = s(tdims(i));
            end
        end
        self.maxNt = max(self.Nts);
        
        % Compute chunksizes for each variable.
        chunksizes = cell(self.Nnames, 1);
        for i = self.allinds
            cs = self.sizes{i};
            if tdims(i) ~= 0
                cs(tdims(i)) = 1;
            end
            chunksizes{i} = cs;
        end
        
        % Do time-varying variables and time-invariant variables.
        m = cell(self.Nnames, 1); % m is short for "map".
        order = {'timevarying', 'timeinvariant'};
        if ~self.parameters_at_end
            order = fliplr(order);
        end
        for o = 1:length(order)
            switch order{o}
            case 'timevarying'
                % Preallocate cell arrays for time-varying variables.
                for i = self.tinds
                    m{i} = cell(self.Nts(i), 1);
                end

                % Loop through time and fields.
                for t = 1:self.maxNt
                    for i = self.tinds
                        if t <= self.Nts(i)
                            [m{i}{t}, v] = self.getchunk(v, chunksizes{i});
                        end
                    end
                end

                % Concatenate chunks.
                for i = self.tinds
                    m{i} = cat(tdims(i), m{i}{:});
                end
            case 'timeinvariant'
                % Variables that don't have a time component.
                for i = self.constinds
                    [m{i}, v] = self.getchunk(v, chunksizes{i});
                end
            otherwise
                error('Invalid operation!');
            end
        end
        
        % Store to the instance.
        self.map = m;
    end%function
    
    function v = vec_numeric(self, s, vtype)
        % Collapses the struct s (in which each entry is a single ND array) into
        % a single vector v using the precomputed indices. This function is
        % used inside a try/catch block, so we don't need to check sizes.
        if isequal(vtype, 'logical') % Fix for Matlab R2015b
            v = false(self.Nvar, 1);
        else
            v = zeros(self.Nvar, 1, vtype);
        end
        for i = self.allinds
            name = self.names{i};
            thismap = self.map{i};
            v(thismap(:)) = s.(name)(:);
        end
    end%function
    
    function v = vec_symbolic(self, s)
        % Collapses the struct s (in which each entry can contain cell arrays)
        % into a single vector v using the precomputed indices.
        vecs = cell(self.Nnames, self.maxNt);
        
        % First do time-varying variables.
        keep = true(size(vecs));
        for i = self.tinds
            name = self.names{i};
            for t = 1:self.Nts(i)
                val = self.indexfinaldim(s.(name), t);
                if iscell(val)
                    val = self.flattencell(val, self.ismat(i));
                elseif self.ismat(i)
                    val = val(:);
                end
                vecs{i, t} = val;
            end
            keep(i, self.Nts(i) + 1:end) = false();
        end
        for i = self.constinds
            keep(i,:) = false();
        end
        vecs = vecs(keep(:));
        
        % Now do constant variables.
        consts = cell(self.Nnames, 1);
        for i = self.constinds
            name = self.names{i};
            consts{i} = s.(name)(:);
        end
        consts = consts(self.constinds);
        
        % Finally, combine everybody into a single vector.
        try
            v = vertcat(vecs{:}, consts{:});
        catch err
            error(['Error concatenating symbolic variables: %s\nDid you ', ...
                   'forget to declare a matrix variable?'], err.message);
        end
    end%function
    
    function [isnum, numtype] = checksize(self, s)
        % Checks that the size of the input struct s matches the original sizes,
        % issuing an error if anything is wrong.
        %
        % Returns a boolean saying whether s is purely numeric or has cell
        % components. Also returns a string of type, defaulting to 'double' if
        % types are mixed and '' if not numeric.
        if ~isstruct(s)
            error('s must be a struct!');
        end
        
        % Make sure fields are the same.
        snames = sort(fieldnames(s));
        if ~isequal(snames, self.names)
            missingnames = setdiff(self.names, snames);
            extranames = setdiff(snames, self.names);
            if ~isempty(missingnames) || ~isempty(extranames)
                error('Invalid fields in s: missing={%s}, extra={%s}', ...
                      mpctools.row2str(missingnames), ...
                      mpctools.row2str(extranames));
            end
        end
        
        % Now check sizes. This one is a bit tricky because we could have some
        % components as cell arrays.
        gotsizes = cell(self.Nnames, 1);
        needsizes = cell(self.Nnames, 1);
        isnum = true();
        numtype = '';
        for i = self.allinds
            name = self.names{i};
            x = s.(name); % Get current variable.
            if self.isarray(x)
                % For arrays, we can actually check the true size.
                gotsizes{i} = size(x);
                needsizes{i} = self.sizes{i};
                if isempty(numtype)
                    numtype = class(x);
                elseif ~isequal(numtype, class(x))
                    numtype = 'double';
                end
            elseif isa(x, 'casadi.SX') || isa(x, 'casadi.MX')
                % Can also just set size, but isn't numeric.
                isnum = false();
                gotsizes{i} = size(x);
                needsizes{i} = self.sizes{i};
            elseif isempty(x) || ~iscell(x)
                error('Invalid entry in field "%s"!', name);
            else
                % For cells, we only compare the non-singleton dimensions. This
                % is not completely general, but it avoids transposition issues.
                isnum = false(); % Set flag.
                
                needsize = self.sizes{i};
                needsize = needsize(needsize ~= 1);
                needsizes{i} = needsize;
                
                gotsize = [size(x{1}), size(x)];
                gotsize = gotsize(gotsize ~= 1);
                gotsizes{i} = gotsize;
            end
        end
        if ~isnum
            numtype = '';
        elseif isempty(numtype)
            numtype = 'double';
        end
        rightsize = cellfun(@isequal, gotsizes, needsizes);
        
        if ~all(rightsize)
            message = cell(self.Nnames, 1);
            for i = self.allinds
                if rightsize(i)
                    message{i} = '';
                else
                    message{i} = sprintf('\n  %s: expect [%s], got [%s]', ...
                                         self.names{i}, ...
                                         mpctools.row2str(needsizes{i}), ...
                                         mpctools.row2str(gotsizes{i}));
                end
            end
            message = horzcat(message{:});
            error('Fields with incorrect sizes: %s', message);
        end
    end%function
    
    function numtype = checknumeric(self, s)
        % Makes sure all fields of s are numeric and returns type. Also checks
        % sizes of each field.
        [isnum, numtype] = self.checksize(s);
        if ~isnum
            % TODO: better error message here.
            error('Some fields of s not numeric! Use symvec for symbols.');
        end
    end%function
end%methods

methods (Static, Access = private)
    function [chunk, v] = getchunk(v, chunksize)
        % Returns the properly-sized chunk from v.
        n = prod(chunksize);
        chunk = reshape(v(1:n), chunksize);
        v = v(n + 1:end);
    end%function

    function cflat = flattencell(c, ismat)
        % Flattens each element of the cell array c and then concatenates.
        if isscalar(c)
            cflat = c{1};
            if ismat
                cflat = cflat(:);
            end
        else
            if ismat
                for i = 1:numel(c)
                    c{i} = c{i}(:);
                end
            end
            cflat = vertcat(c{:});
        end
    end%function
    
    function xt = indexfinaldim(x, t)
        % Returns `x(...,t)`, where `...` is an appropriate number of colons.
        switch ndims(x)
        case 2
            xt = x(:,t);
        case 3
            xt = x(:,:,t);
        otherwise
            idx = substruct('()', [repmat({':'}, 1, ndims(x) - 1), {t}]);
            xt = subsref(x, idx);     
        end
    end%function
    
    function tf = isarray(x)
        % Returns true or false whether x is an array of doubles, logicals, or
        % characters.
        tf = isnumeric(x) || islogical(x) || ischar(x);
    end%function
end%methods

methods (Static, Access = public)
    function names = sortnames(names)
        % Sorts variables in an appropriate order for MPC.
        frontvars = {'x'; 'xc'; 'y'; 'z'; 'zc'; 'v'; 'u'; 'w'};
        [found, foundind] = ismember(names, frontvars);
        foundind = foundind(foundind ~= 0);
        names = [frontvars(sort(foundind)); sort(names(~found))];
    end%function
end%methods

end%classdef

