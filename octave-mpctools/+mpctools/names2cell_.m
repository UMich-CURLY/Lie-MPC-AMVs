function names = names2cell_(namesmat)
% names = names2cell_(namesmat)
%
% Converts a matrix of variable names to a row cell array.
%
% namesmat should be a char matrix, e.g., the output of Casadi's
% Function.name_in() or Function.name_out().
narginchk(1, 1)
namesmat(namesmat == 0) = ' '; % Required for Matlab.
nrows = size(namesmat, 1);
ncols = size(namesmat, 2);
names = strtrim(mat2cell(namesmat, ones(nrows, 1), ncols)');
end%function

