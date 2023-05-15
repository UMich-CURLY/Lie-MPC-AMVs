function [Ainv, Achol] = spdinv(A)
% `[Ainv, Achol] = spdinv(A)`
%
% Computes inverse of symmetric positive-definite matrix `A` via
% (upper-triangular) Cholesky factorization. `A` must be given as a single
% positional argument.
%
% Also returns the Cholesky factor `Achol` (with `Achol'*Achol == A`).
%
% Note that `A` can be a 3D array, at which point a 3D array of inverses is
% returned assuming each each `A(:,:,i)` slice is a matrix.
I = eye(size(A, 1));
switch ndims(A)
case 2
    [Ainv, Achol] = spdinv_(A, I);
case 3
    Ainv = NaN(size(A));
    Achol = NaN(size(A));
    for i = 1:size(A, 3)
        [Ainv(:,:,i), Achol(:,:,i)] = spdinv_(A(:,:,i), I);
    end
otherwise
    error('A must be a square matrix or a 3D array!');
end
end%function

function [Ainv, Achol] = spdinv_(A, I)
    % Calculates inverse for a single matrix.
    Achol = chol(A);
    Ainv = Achol\(Achol'\I);
end%function
