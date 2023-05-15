classdef TestScript < handle
% Class for running and holding results of a test script.

properties (GetAccess=public, SetAccess=private)
    script;
    fails = {};
    successes = {};
    runtimeerror = [];
    output;
    workspace;
end

methods (Access=public)
    function self = TestScript(scriptname)
        % TestScript(scriptname)
        %
        % Runs the given script, changing any assert() calls to Assertion. If
        % the file fails for some other reason, that error is stored as
        % runtimeerror.
        narginchk(1, 1);
        self.script = scriptname;
        self.workspace = self.execscript();
    end%function
    
    function print(self, printmode)
        % self.print(['quiet'|'normal'|'verbose'])
        %
        % Prints the results of the script.
        narginchk(1, 2);
        if nargin() < 2
            printmode = 'normal';
        end
        printsuccess = true();
        printverbose = false();
        switch printmode
        case 'quiet'
            printsuccess = false();
        case 'normal'
            % Pass.
        case 'verbose'
            printverbose = true();
        otherwise
            error('Invalid argument!');
        end
        
        % Stop early if in quiet mode.
        if ~printsuccess && self.issuccessful()
            return
        end
        
        % Print results.
        fprintf('*** Test results for script "%s":\n', self.script);
        if printsuccess
            fprintf('    %d Successes\n', length(self.successes));
        end
        fprintf('    %d Failures\n', length(self.fails));
        if printverbose
            for i = 1:length(self.fails)
                a = self.fails{i};
                fprintf('*** Assertion failed on line %d:\n', ...
                        a.stack(end).line);
                self.printerrormessage(a.message);
            end
            fprintf('*** Script output:\n');
            self.printerrormessage(self.output);
        end
        if ~isempty(self.runtimeerror)
            fprintf('!!! Uncaught error on line %d:\n', ...
                    self.runtimeerror.stack(end).line);
            self.printerrormessage(self.runtimeerror.message);
        end
        if self.issuccessful()
            status = 'Success';
        else
            status = 'FAILURE';
        end
        fprintf('*** %s\n\n', status);
    end%function
    
    function tf = issuccessful(self)
        % Returns whether test was successful.
        tf = isempty(self.runtimeerror) && isempty(self.fails);
    end%function
end%methods

methods (Access=private)
    function ws = execscript(self)
        % Executes the current script, populating assertions and errors.
        assert = @(varargin) self.assertion(varargin{:});
        try
            self.output = evalc(self.script);
        catch err
            self.runtimeerror = err;
        end
        workspacevars = setdiff(who(), {'assert', 'self'});
        ws = struct();
        for i = 1:length(workspacevars)
            v = workspacevars{i};
            ws.(v) = eval(v);
        end
    end%function
    
    function assertion(self, varargin)
        % Stores assertion to the current object.
        a = Assertion(varargin{:});
        a.stack = a.stack(1:(end - 3));
        if a.successful
            self.successes{end + 1} = a;
        else
            self.fails{end + 1} = a;
        end
    end%function
end%methods

methods (Static, Access=private)
    function printerrormessage(message)
        narginchk(1, 1);
        messagelines = strsplit(message, sprintf('\n'));
        if length(messagelines) > 1 && isempty(messagelines{end})
            messagelines = messagelines(1:(end - 1));
        end
        fprintf('    | %s\n', messagelines{:});
    end%function
end%methods

end%classdef

