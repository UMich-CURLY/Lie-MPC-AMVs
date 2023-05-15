% Checks to make sure Matlab implementation of colloc matches Octave's
% (correct) implementation.
assert(mpctools.isOctave()); % This test only runs on Octave.
Ns = 2:50;
errs = NaN(length(Ns), 4);
fprintf('log10 of maximum error (should all be small)\n');
fprintf('%2s : %10s %10s %10s %10s\n', 'N', 'r', 'A', 'B', 'q');
fprintf('%s\n', repmat('-', 1, 75));
for i = 1:length(Ns)
    N = Ns(i);
    calc = cell(4, 1);
    [calc{:}] = mpctools.collocweights(N, 'matlab');
    
    check = cell(4, 1);
    [check{:}] = colloc(N);
    
    errs(i,:) = cellfun(@(a, b) log10(max(abs(a(:) - b(:)))), calc, check);
    fprintf('%2d : %10g %10g %10g %10g\n', N, errs(i,:));
end

% Assert bounds.
assert(all(errs(:,1) < -12)); % Error on roots.
assert(all(errs(:,2) < -8)); % Error on 1st derivative weights.
assert(all(errs(:,3) < -5)); % Error on 2nd derivative weights.
assert(all(errs(:,4) < -10)); % Error on quadrature weights.

