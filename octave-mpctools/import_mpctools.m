function varargout = import_mpctools(local)
% `mpctools = import_mpctools()`
% `import_mpctools('*')`
%
% Returns a "module" struct with mpctools functions. This is similar to the
% Python construct `import mpctools`
%
% If the second form is used, the function handles are also created in the
% calling namespace. This is similar to the Python construct
% `from mpctools import *`. Note that we do not advise using this form, as it
% is only for lazy people who like namespace pollution.
narginchk(0, 1);
if nargin() > 0 && isequal(local, '*')
    importstar = true();
    nargoutchk(0, 0);
elseif nargin() == 0
    importstar = false();
    nargoutchk(0, 1);
else
    error('Unknown argument: %s.', disp(local));
end
importlist = { ...
    'c2d';
    'collocweights';
    'ekf';
    'getCasadiFunc';
    'getCasadiIntegrator';
    'getCasadiDAE';
    'getLinearizedModel';
    'KalmanFilter';
    'rk4';
    'mpcplot';
    'nlfilter';
    'nmpc';
    'nmhe';
    'smushcolloc';
    'spdinv';
    'sstarg';
    'jacobian';
    'version';
};

module = struct();
for i = 1:length(importlist)
    func = importlist{i};
    try
        module.(func) = str2func(['mpctools.', func]);
    catch err
        error('Error importing %s: %s', func, err.message);
    end
    if importstar
        assignin('caller', func, module.(func));
    end
end

if importstar
    varargout = {};
else
    varargout = {module};
end

end%function

