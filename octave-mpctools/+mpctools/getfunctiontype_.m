function casaditype = getfunctiontype_(func)
% casaditype = getfunctiontype_(func)
%
% Returns the CasADi type ('SX' or 'MX') of the given function object.
narginchk(1, 1);
if ~mpctools.iscasadifunc(func)
    error('func is not a casadi.Function object!');
end
typename = func.type_name();
if isequal(typename, 'Function')
    typename = func.class_name();
end
switch typename
case {'sxfunction', 'SXFunction'}
    casaditype = 'SX';
case {'mxfunction', 'MXFunction'}
    casaditype = 'MX';
otherwise
    error('Unknown function type "%s"!', typename);
end
end%function

