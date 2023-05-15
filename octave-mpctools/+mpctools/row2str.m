function s = row2str(row, fmt, sep)
% s = row2str(row, [fmt], [sep])
%
% Returns a row vector as a one-line string (without brackets or braces).
% Second argument fmt specifies the printf formatting string. By default, it is
% guessed from the input. Third argument sep is the separator string.
% the default is ', '.
narginchk(1, 3);
if nargin() < 2
    if isnumeric(row)
        fmt = '%g';
    else
        fmt = '%s';
    end
end
if nargin() < 3
    sep = ', ';
end
row = row(:); % Flatten any arrays.
fmt = repmat([fmt, sep], 1, length(row));
if ~isempty(fmt)
    fmt = fmt(1:(end - length(sep)));
end
if iscell(row)
    s = sprintf(fmt, row{:});
else
    s = sprintf(fmt, row);
end

end%function

