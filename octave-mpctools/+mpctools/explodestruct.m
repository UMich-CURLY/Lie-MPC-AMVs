function c = explodestruct(s)
% c = explodestruct(s)
%
% Explodes struct s into cell array of {field, value} pairs.
narginchk(1, 1);
if ~isstruct(s)
    error('Input s must be a struct!');
end
fields = fieldnames(s);
c = cell(2*length(fields), 1);
for i = 1:length(fields)
    f = fields{i};
    c{2*i - 1} = f;
    c{2*i} = s.(f);
end
end%function

