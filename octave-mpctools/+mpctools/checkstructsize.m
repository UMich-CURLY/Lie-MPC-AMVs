function checkstructsize(struct_, sizes_)
% checkstructsize(struct, sizes)
%
% Checks whether structsize(struct) == sizes. If not, issues an error.
missingfields = setdiff(fieldnames(sizes_), fieldnames(struct_));
if ~isempty(missingfields)
    fmt = repmat(' "%s"', 1, length(missingfields));
    error(['Missing fields in %s:', fmt], inputname(1), missingfields{:});
end
checksize = mpctools.structsize(struct_);
fieldsokay = struct();
fields = fieldnames(struct_);
badfields = cell(size(fields));
for i = 1:length(fields)
    f = fields{i};
    fieldsokay.(f) = isequal(checksize.(f), sizes_.(f));
    if ~fieldsokay.(f)
        badfields{i} = sprintf('"%s" (is [%s], should be [%s])', f, ...
                               mpctools.row2str(checksize.(f)), ...
                               mpctools.row2str(sizes_.(f)));
    end
end
fieldsokay = mpctools.explodestruct(fieldsokay);
fields = fieldsokay(1:2:end);
okay = cell2mat(fieldsokay(2:2:end));
if ~all(okay)
    badfields = badfields(~okay);
    error('Incorrect sizes in %s: %s', inputname(1), ...
          mpctools.row2str(badfields, '%s', '; '));
end
end%function

