classdef ArgumentParser < handle
% Parser for keyword-style arguments.

properties (Access=private)
    argnames_ = {};
    checkers_ = struct();
    defaults_ = struct();
    required_ = {};
    DEFAULT_VALUE = [];
    ALLOW_UNKNOWN_KWARGS = false();
    ALLOW_DUPLICATE_KW = true();
end%properties

methods (Access=public)
    function add(self, argname, default, checker)
        % self.add(argname, [default], [checker])
        %
        % Adds an argument, optinally with a default value.
        %
        % Second argument checker is a function that takes one argument and
        % returns a (possibly modified) version, raising an error if the
        % input is invalid. Note that these can also be any strings accepted
        % by ArgumentParser.checkarg().
        narginchk(2, 4);
        if ~ischar(argname)
            error('argname must be a string!');
        elseif isequal(argname, '**kwargs')
            self.kwargs(true()); % Allow keyword arguments.
            return
        elseif ~isvarname(argname) && ~iskeyword(argname)
            error('Invalid argument name: "%s"', argname);
        elseif ismember(argname, self.argnames_)
            error('Duplicate argument name: "%s"', argname);
        end
        self.argnames_ = [self.argnames_; {argname}];
        if nargin() < 3 || isequal(default, 'required')
            self.required_ = [self.required_; {argname}];
        else
            self.defaults_.(argname) = default;
        end
        if nargin() >= 4 && ~isempty(checker)
            self.checkers_.(argname) = checker;
        end
    end%function
    
    function required = get_required(self)
        % self.get_required()
        %
        % Returns a cell array of required arguments.
        required = self.required_;
    end%function
    
    function defaults = get_defaults(self)
        % self.get_defaults()
        %
        % Returns struct of default values.
        defaults = self.defaults_;
    end%function
    
    function [args, extra, isdefault_] = parse(self, varargin)
        % [args, extra, isdefault] = self.parse(...)
        %
        % parses arguments based on given values.
        %
        % Arguments can be positional, (key, value) pairs. The key '-struct' or
        % '**' indicates that the corresponding value is a struct whose fields
        % should be treated as keyword arguments.
        %
        % Second output extra is a struct with all the unknown arguments. It is
        % empty unless '**kwargs' has been added as an argument.
        %
        % Third output isdefault is a struct that specifies whether each
        % argument is the default value or not.
        args = struct();
        extra = struct();
        N = length(varargin);
        Npos = min(N, length(self.argnames_));
        n = 1;
        
        % Add positional arguments.
        while n <= Npos && ~ischar(varargin{n})
            args.(self.argnames_{n}) = varargin{n};
            n = n + 1;
        end
        
        % Add keyword arguments.
        while n <= N
            kw = varargin{n};
            if ~ischar(kw) || n == N
                error('keyword arguments must come in pairs!');
            elseif ismember(kw, {'-struct', '**'})
                % Struct. Process all fields.
                thisstruct = varargin{n + 1};
                if ~isstruct(thisstruct)
                    error('Value following "**" must be a struct!');
                end
                structkw = fieldnames(varargin{n + 1});
                for i = 1:length(structkw)
                    kw = structkw{i};
                    [args, extra] = self.parsekw(kw, thisstruct.(kw), ...
                                                 args, extra);
                end
            else
                % Normal keyword argument.
                [args, extra] = self.parsekw(kw, varargin{n + 1}, args, extra);
            end
            n = n + 2;
        end
        
        % Add in defaults and check given values.
        isdefault_ = struct();
        for i = 1:length(self.argnames_)
            kw = self.argnames_{i};
            isdefault_.(kw) = ~isfield(args, kw) || self.isdefault(args.(kw));
            if ~isdefault_.(kw)
                % Check value.
                if isfield(self.checkers_, kw)
                    checker = self.checkers_.(kw);
                    args.(kw) = self.checkarg(kw, args.(kw), checker);
                end
            elseif isfield(self.defaults_, kw)
                % Use default value.
                args.(kw) = self.defaults_.(kw);
            elseif ismember(kw, self.required_)
                % Required argument not given.
                error('Argument "%s" is required!', kw);
            else
                % Fall back to default default.
                args.(kw) = self.DEFAULT_VALUE;
            end
        end
    end%function
    
    function kwargs(self, tf)
        % self.kwargs(true())
        % self.kwargs(false())
        %
        % Allow or dis-allow unknown keyword arguments. If allowed, they will be
        % returned in the extra argument of self.parse().
        %
        % Note that self.add('**kwargs') is equivalent to self.kwargs(true()).
        if ~isscalar(tf) || ~islogical(tf)
            error('Input must be true or false!');
        end
        self.ALLOW_UNKNOWN_KWARGS = tf;
    end%function
end%methods

methods (Access=private)
        function [args, extra] = parsekw(self, kw, val, args, extra)
        % self.parsekw(kw, val, args, extra)
        %
        % Adds a keyword argument to args or (depending on whether unknown
        % arguments are allowed) adds it to extra/raises an error.
        narginchk(5, 5);
        if ismember(kw, self.argnames_)
            if ~self.ALLOW_DUPLICATE_KW
                if isfield(args, kw)
                    error('Duplicate keyword "%s".', kw);
                end
            end
            args.(kw) = val;
        elseif self.ALLOW_UNKNOWN_KWARGS
            extra.(kw) = val;
        else
            error('Unknown keyword argument: "%s"', kw);
        end
    end%function
    
    function tf = isdefault(self, val)
        % tf = self.isdefault(val)
        %
        % Returns true or false whether the value is exactly equal to the
        % default sentinal value.
        tf = isequal(val, self.DEFAULT_VALUE) ...
             && isequal(class(val), class(self.DEFAULT_VALUE));
    end%function
end%methods

methods (Static)
    function x = checkarg(name, x, checks)
        % x = checkarg(name, x, checks)
        %
        % Checks an argument and possibly fixes it. checks can be a string from
        % the following list:
        %
        % - 'row', 'col' (allows any vector but reshapes to row or col)
        % - 'matrix'
        % - 'cell', 'numeric', 'struct'
        % - 'scalar'
        % - 'pos', 'nonneg'
        % - 'int', 'logical'
        % - 'string' (allows a string or 1-element cell array with a string)
        % - 'varname'
        % - 'bool' (scalar logical)
        % - 'casadifunc' (casadi.Function object)
        % - 'func' (function handle or casadi.Function)
        % - '<...>' (object of type ...)
        % - '/.../' (string exactly matching regex ...)
        %
        % checks can also be a cell array of strings to apply all of the given
        % checks (in the given order)
        narginchk(3, 3);
        if ~iscell(checks)
            checks = {checks};
        end
        % Build array with simple "is*" checks.
        is = struct('cell', @iscell, 'numeric', @isnumeric, ...
                    'struct', @isstruct, 'varname', @isvarname, ...
                    'scalar', @isscalar, 'matrix', @ismatrix, ...
                    'logical', @islogical);
        % Now apply all the checks.
        for i = 1:length(checks)
            check = checks{i};
            if isa(check, 'function_handle')
                try
                    x = check(x);
                catch err
                    error('Argument %s is invalid: %s', name, err.message);
                end
                continue
            end
            switch check
            case 'row'
                if ~isvector(x) && ~isempty(x)
                    error('Argument %s must be a vector!', name);
                end
                x = reshape(x, [1, length(x)]);
            case {'col', 'column'}
                if ~isvector(x) && ~isempty(x)
                    error('Argument %s must be a vector!', name);
                end
                x = reshape(x, [length(x), 1]);
            case {'pos', 'positive'}
                if ~all(x(:) > 0)
                    error('Argument %s must be strictly positive!', name);
                end
            case {'nonneg', 'nonnegative'}
                if ~all(x(:) >= 0)
                    error('Argument %s must be nonnegative!', name);
                end
            case {'int', 'integer'}
                if ~all(round(x(:)) == x(:))
                    error('Argument %s must have integer values!', name);
                end
            case {'str', 'string'}
                % String arguments can be passed as one-element cell arrays to
                % avoid being interpreted as a keyword argument.
                if iscell(x)
                    x = x{1};
                end
                if ~ischar(x)
                    error('Argument %s must be a string!', name);
                end
            case {'bool', 'boolean'}
                if ~isscalar(x) || ~islogical(x)
                    error('Argument %s must be a (scalar) boolean!', name);
                end
            case 'casadifunc'
                if ~mpctools.iscasadifunc(x)
                    error('%s must be a casadi.Function!', name);
                end
            case 'func'
                if ~isa(x, 'function_handle') && ~mpctools.iscasadifunc(x)
                    error(['%s must be a function handle ', ...
                           '(or a casadi.Function)!'], name);
                end
            otherwise
                if isfield(is, check)
                    % Simple checks.
                    if ~is.(check)(x)
                        error('Argument %s must be %s!', name, check);
                    end
                elseif mpctools.ArgumentParser.isdelimited(check, '<>')
                    % Custom classes.
                    if ~isa(x, check(2:(end - 1)))
                        error('Argument %s must be a %s object!', name, check);
                    end
                elseif mpctools.ArgumentParser.isdelimited(check, '//')
                    % Regex match.
                    x = mpctools.ArgumentParser.checkarg(name, x, 'string');
                    pattern = check(2:(end - 1));
                    if isempty(regexp(x, ['^', pattern, '$'], 'once'))
                        error('Argument %s must match "%s"!', name, pattern);
                    end
                else
                    error('Invalid check for argument %s: "%s"', name, check)
                end
            end
        end 
    end%function
    
    function [tf, contents] = isdelimited(s, d1, d2)
        % isdelimited(s, delimiters)
        %
        % Checks whether the string s is delimited by d1 and d2. If d1 and d2
        % are each single characters, they can be passed as a single string,
        % e.g., '<>' for angle brackets or '()' for parentheses.
        narginchk(2, 3);
        if nargin() < 3
            if length(d1) == 2
                d2 = d1(2);
                d1 = d1(1);
            else
                error('Invalid delimiter specified!');
            end
        end
        tf = length(s) >= length(d1) + length(d2) ...
             && strcmp(d1, s(1:length(d1))) ...
             && strcmp(d2, s((end - length(d2) + 1):end));
         if tf
            contents = s((length(d1) + 1):(end - length(d2)));
         else
            contents = '';
         end
    end%function
end%methods

end%classdef

