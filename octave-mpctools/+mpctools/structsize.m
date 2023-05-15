function ss = structsize(s)
% ss = structsize(s)
%
% Returns a struct ss whose elements are the sizes of the elements of s.
ss = structfun(@size, s, 'UniformOutput', false());
end%function

