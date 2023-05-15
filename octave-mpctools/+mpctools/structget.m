function [d, s, ffound] = structget(s, f, d)
% [val, ss, ffound] = structget(s, f, d)
%
% returns s.(f) if f is a field name for s. Otherwise, returns d. Also returns
% ss, which is a copy of s with field f removed (if present). Also returns a
% boolean value ffound to say whether the field f was found in struct s.

ffound = isfield(s, f);
if ffound
    d = s.(f);
    if nargout() > 1
        s = rmfield(s, f);
    end
end

end%function
