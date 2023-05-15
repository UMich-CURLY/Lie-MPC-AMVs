function fnew = changefuncargs_(f, newargs)
% fnew = changefuncargs_(f, newargs)
%
% Returns a new casadi.Function object whos inputs have the given argument
% names.
narginchk(2, 2);
if ~mpctools.iscasadifunc(f)
    error('f must be a casadi.Function!');
elseif length(newargs) ~= f.n_in()
    error('f expects %d inputs, but %d were given!', f.n_in(), length(newargs));
end
sizes = cell(size(newargs));
for i = 1:length(sizes)
    sizes{i} = f.size_in(i - 1);
end
fnew = mpctools.getCasadiFunc(f, sizes, newargs, {f.name_out(0)}, ...
                              'casaditype', mpctools.getfunctiontype_(f));
end%function

