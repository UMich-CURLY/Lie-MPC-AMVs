classdef Assertion < handle
% Runs an assert statement but doesn't error on failure.

properties (Access=public)
    stack;
    args;
    successful;
    message;
end%properties

methods (Access=public)
    function self = Assertion(varargin)
        % Assertion(...)
        %
        % Initialize the object as you would an assert statement. The condition
        % will be checked, but no error will be issued on failure.
        self.stack = self.getstack();
        self.args = varargin;
        try
            assert(varargin{:});
            self.message = '';
            self.successful = true();
        catch err
            self.message = err.message;
            self.successful = false();
        end
    end%function
end%methods

methods (Static)
    function s = getstack(n)
        % s = UnitTest.getstack([n=0])
        %
        % Returns current stack trace, removing the final n frames.
        narginchk(0, 1);
        if nargin() < 2
            n = 0;
        end
        try
            error(' ');
        catch err
            s = err.stack;
        end
        s = s(1:(end - n));
    end%function
end%methods

end%classdef

