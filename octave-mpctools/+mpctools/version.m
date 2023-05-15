function varargout = version()
% `version()`
%
% Returns the version of MPCTools as a string.
%
% The second output argument returns the HG changeset ID (as a hexidecimal
% string), which is more granular than the version number.
narginchk(0, 0);
versionstring = '2.3.3';
hgid = '%HG_CHANGESET_ID%';
% <--
thisdir = fileparts(mfilename('fullpath'));
mpctoolsdir = fileparts(thisdir);
[status, hgid] = system(sprintf('hg id -i ''%s''', mpctoolsdir));
if status ~= 0
    warning('Error running `hg id` (exited with status %d)!', status);
end
hgid = strtrim(hgid);
% -->
if nargout() == 0
    fprintf('MPCTools Version %s (hg changeset %s)\n', versionstring, hgid);
    varargout = {};
else
    varargout = {versionstring, hgid};
end
end%function

