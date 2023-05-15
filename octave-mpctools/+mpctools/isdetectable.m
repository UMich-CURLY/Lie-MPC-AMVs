function [tf, Nundetect, Abar, Bbar, Cbar, T] = isdetectable(varargin)
% `[bool, Nundetect, Abar, Bbar, Cbar, T] = isdetectable(A, B, C, [tol])`
%
% Determines whether a system is detectable or not.
%
% The system is first put in observability canonical form, and then the
% eigenvalues of the unobservable modes are checked.
%
% If provided, `tol` is a scalar tolerance for stability. Unobservable
% eigenvalues must have absolute value below `1 - tol` for the system to be
% considered detectable.
persistent parser
if isempty(parser)
    parser = mpctools.ArgumentParser();
    parser.add('A', 'required', 'numeric');
    parser.add('B', 'required', 'numeric');
    parser.add('C', 'required', 'numeric');
    parser.add('tol', 0, {'scalar', 'numeric'});
end
args = parser.parse(varargin{:});

% Special case in Octave.
if nargout() == 1 && mpctools.isOctave()
    tf = builtin('isdetectable', args.A, args.B, args.C, args.tol);
    return
end

% Method that works for Matlab.
[Abar, Bbar, Cbar, T, Ndetect] = obsvf(args.A, args.B, args.C);
Ndetect = sum(Ndetect);

Ano = Abar((Ndetect + 1):end,(Ndetect + 1):end);
noeigs = eig(Ano);
Nundetect = nnz(abs(noeigs) >= 1 - args.tol);
tf = (Nundetect == 0);

end%function

