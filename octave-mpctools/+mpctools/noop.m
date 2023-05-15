function varargout = noop(varargin)
% varargout = noop(varargin)
%
% Function that does nothing. Can be used as a callback placeholder or to
% eliminate other functions, e.g. `fprintf = @noop` to shut off printing.
%
% All return values are 0 by 0 empty matrices.
varargout = repmat({[]}, 1, nargout());
end%function
