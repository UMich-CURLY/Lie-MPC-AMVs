function verb = MAX_VERBOSITY(val)
% verb = MAX_VERBOSITY([val])
%
% Returns the current state of the maximum verbosity and optionally sets it
% to a new value.
persistent maxverbosity
if isempty(maxverbosity)
    maxverbosity = 100;
end
narginchk(0, 1);
verb = maxverbosity;
if nargin() > 0
    maxverbosity = val;
end
end%function

