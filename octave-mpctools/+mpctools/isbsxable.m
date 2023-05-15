function [tf, r1, r2] = isbsxable(s1, s2)
% [tf, r1, r2] = isbsxable(s1, s2)
%
% Returns true or false if the two shape vectors are compatable for binary
% singleton expansion. Note that they are 1-padded to be the same length.
%
% Additional outputs r1 and r2 are the "repmat sizes" for the two inputs that
% would give the full-size array if passed to repmat. If tf is false, both
% of these are NaNs.

nd = max(length(s1), length(s2));
s1 = [s1, ones(1, nd - length(s1))];
s2 = [s2, ones(1, nd - length(s2))];

s1one = (s1 == 1);
s2one = (s2 == 1);
tf = all(s1one | s2one | (s1 == s2));

if tf
    r1 = ones(size(s1));
    r1(s1one) = s2(s1one);
    
    r2 = ones(size(s2));
    r2(s2one) = s1(s2one);
else
    r1 = NaN(size(s1));
    r2 = NaN(size(s2));
end

end%function
