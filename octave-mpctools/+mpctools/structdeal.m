function varargout = structdeal(s, varargin)
% [val1, val2, ...] = structdeal(struct, field1, field2, ...)
% [val1, val2, ...] = structdeal({struct, default}, field1, field2, ...)
%
% Deals the given fields of the struct.
%
% By default an error is issued if any field are not found. However, is s is
% given as a 2-element cell array, then the default value is used for any
% fields that don't exist.
if iscell(s)
    if numel(s) ~= 2
        error('Invalid input for s!');
    end
    usedefault = true();
    default = s{2};
    s = s{1};
else
    usedefault = false();
end
if ~isstruct(s)
    error('s must be a struct!');
end
varargout = cell(size(varargin));
for i = 1:length(varargin)
    f = varargin{i};
    if ~ischar(f)
        error('field%s is not a string!', i);
    end
    if isfield(s, f)
        varargout{i} = s.(f);
    elseif usedefault
        varargout{i} = default;
    else
        error('Invalid output field "%s"!', f);
    end
end
end%function
