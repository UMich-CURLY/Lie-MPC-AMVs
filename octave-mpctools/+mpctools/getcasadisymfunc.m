function func = getcasadisymfunc(casaditype)
% func = getcasadisymfunc(casaditype)
%
% Returns a function handle to casadi.SX.sym or casadi.MX.sym depending on
% whether casaditype is 'SX' or 'MX'.
%
% Note that func is an anonymous function handle due to Octave's lack of
% support for direct function handles to package or static class functions.
if ~ischar(casaditype)
    error('casaditype must be a string!');
end
switch upper(casaditype)
case 'SX'
    func = @(varargin) casadi.SX.sym(varargin{:});
case 'MX'
    func = @(varargin) casadi.MX.sym(varargin{:});
otherwise
    error('Unknown value for casaditype: %s', casaditype);
end

end%function

