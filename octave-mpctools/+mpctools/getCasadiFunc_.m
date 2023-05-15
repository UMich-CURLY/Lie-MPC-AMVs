function things = getCasadiFunc_(f, varsizes, varnames, funcname, casaditype)
% things = getCasadiFunc_(f, varsizes, varnames, funcname, casaditype)
%
% Logic behind getCasadiFunc. Returns a struct with fields fexpr, args, names,
% and sizes.
narginchk(5, 5);

% Figure out sizes and names.
if isvector(varsizes) && ~iscell(varsizes)
    varsizes = num2cell(varsizes);
elseif ~iscell(varsizes)
    error('Invalid input for varsizes!');
end
nvar = length(varsizes);
if isempty(varnames)
    varnames = cell(nvar, 1);
    for i = 1:nvar
        varnames{i} = sprintf('x_%d', i);
    end
elseif ~iscell(varnames)
    error('Invalid input for varnames!');
elseif length(varsizes) ~= length(varnames)
    error('Sizes and names must be the same length!');
end
if ~ischar(funcname) || ~isvarname(funcname)
    error('Invalid input for funcname!');
end

% Make sure cells are rows.
varnames = varnames(:)';
varsizes = varsizes(:)';

% Choose symbolic function.
symfunc = mpctools.getcasadisymfunc(casaditype);

% Create variables and function expression.
inputvar = cell(1, nvar); % Needs to be row cell array.
for i = 1:nvar
    inputvar{i} = symfunc(varnames{i}, round(varsizes{i}));
end
fexpr = f(inputvar{:});

% Build return struct.
things = struct();
things.fexpr = fexpr;
things.args = inputvar;
things.names = varnames;
things.sizes = varsizes;

end%function
