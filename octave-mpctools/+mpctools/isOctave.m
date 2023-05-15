function octavetf = isOctave(check)
% octavetf = isOctave()
%
% Checks to see if script/function is being called by Octave.
%
% This is done by looking for a builtin function called "OCTAVE_VERSION", so
% if you somehow add this function to Matlab, then there will be problems.

persistent OCTAVETF

if isempty(OCTAVETF) || nargin() > 0
	OCTAVETF = logical(exist('OCTAVE_VERSION','builtin'));
end

octavetf = OCTAVETF;

end
