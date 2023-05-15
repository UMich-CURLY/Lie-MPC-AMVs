function varargout = rununittests(runmode)
% rununittests(['quiet'|'normal'|'verbose'])
% tests = rununittests(...)
% [failures, successes] = rununittests(...)
%
% Runs all test scripts in the current directory and returns a struct of
% TestScript objects.
%
% If no outputs are requested, function will error at the end if any tests
% fail. Otherwise, test results will be returned as structs.
narginchk(0, 1);
nargoutchk(0, 2);

% Check input argument.
if nargin() < 1
    runmode = 'normal';
end
quiet = isequal(runmode, 'quiet');

% Run each test script.
files = { ...
    'testc2d';
    'testcolloc';
    'testdae';
    'testgetfunctiontype';
    'testgetqp';
    'testlinearizedmodel';
    'testrk4';
    'testsingleshooting';
    'teststructvect';
};
tests = struct();
for i = 1:length(files)
    f = files{i};
    if ~quiet
        fprintf('Running %s\n', f);
    end
    tests.(f) = TestScript(f);
end

% Print status for each.
if ~quiet
    fprintf('\n');
end
failures = struct();
successes = struct();
for i = 1:length(files)
    f = files{i};
    t = tests.(f);
    t.print(runmode);
    if t.issuccessful()
        successes.(f) = t;
    else
        failures.(f) = t;
    end
end

% Print final status method.
Nfailures = length(fieldnames(failures));
if Nfailures == 0
    if ~quiet
        fprintf('All tests successful.\n');
    end
elseif nargout() > 0
    fprintf('%d tests failed!\n', Nfailures);
end

% Choose output arguments.
switch nargout()
case 0
    varargout = {};
    if Nfailures > 0
        error('%d tests failed!', Nfailures);
    end
case 1
    varargout = {tests};
case 2
    varargout = {successes, failures};
end

end%function

