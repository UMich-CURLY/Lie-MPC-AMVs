function s = structupdate(s, u, keys)
% s = structupdate(s, u, [keys])
%
% Updates the fields of struct s using the fields from struct u. If cell array
% keys is provided, only update the fields given in keys.
narginchk(2, 3);
if nargin() < 3
    keys = fieldnames(u);
end
for i = 1:length(keys)
    f = keys{i};
    s.(f) = u.(f);
end

end%function
