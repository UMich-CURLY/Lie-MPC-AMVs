% Unit tests for c2d.
mpc = import_mpctools();

function mats = oldc2d(args)
    % Old method from Pannocchia, Rawlings, Mayne, and Mancuso (2014). Not
    % great if A or -A has large eigenvalues.
    [Nx, Nu] = size(args.b);
    
    % Build weird large matrix thing.
    i = cell(4, 1);
    for j = 1:3
        i{j} = ((j - 1)*Nx + 1):(j*Nx);
    end
    i{4} = (3*Nx + 1):(3*Nx + Nu);
    M = zeros(3*Nx + Nu, 3*Nx + Nu);
    M(i{1},i{1}) = -args.a';
    M(i{2},i{2}) = -args.a';
    M(i{3},i{3}) = args.a;
    M(i{1},i{2}) = eye(Nx);
    M(i{2},i{3}) = args.q;
    M(i{3},i{4}) = args.b;
    
    % Exponentiate and split into parts.
    M = expm(args.Delta*M);
    F3 = M(i{3},i{3});
    G3 = M(i{3},i{4});
    G2 = M(i{2},i{3});
    H2 = M(i{2},i{4});
    K1 = M(i{1},i{4});
    
    % Use formulas.
    mats.A = F3;
    mats.B = G3;
    mats.Q = F3'*G2;
    mats.M = F3'*H2;
    mats.R = args.r*args.Delta + args.b'*F3'*K1 + K1'*F3*args.b;
end%function

function mats = lessoldc2d(args)
    % Alternative method from "Computing Integrals involving the Matrix
    % Exponential" (Van Loan, 1978). Doesn't include a continuous-time m
    % term.
    [Nx, Nu] = size(args.b);
    
    % Build weird large matrix thing.
    i = cell(4, 1);
    for j = 1:3
        i{j} = ((j - 1)*Nx + 1):(j*Nx);
    end
    i{4} = (3*Nx + 1):(3*Nx + Nu);
    M = zeros(3*Nx + Nu, 3*Nx + Nu);
    M(i{1},i{1}) = -args.a';
    M(i{2},i{2}) = -args.a';
    M(i{3},i{3}) = args.a;
    M(i{1},i{2}) = eye(Nx);
    M(i{2},i{3}) = args.q;
    M(i{3},i{4}) = args.b;
    
    % Exponentiate and grab pieces.
    [M, nsquare] = expmnosquare(M*args.Delta);
    F3 = M(i{3},i{3});
    G3 = M(i{3},i{4});
    G2 = M(i{2},i{3});
    H2 = M(i{2},i{4});
    K1 = M(i{1},i{4});
    
    % Do recursive formulas for repeated squaring.
    A = F3;
    B = G3;
    Q = F3'*G2;
    M = F3'*H2;
    R = args.b'*F3'*K1 + K1'*F3*args.b;
    for k = 1:nsquare
        R = 2*R + B'*M + M'*B + B'*Q*B;
        M = M + A'*(Q*B + M);
        Q = Q + A'*Q*A;
        B = B + A*B;
        A = A^2;
    end
    
    % Package up matrices, adding extra part for R.
    mats = struct('A', A, 'B', B, 'Q', Q, 'M', M, 'R', R + args.r*args.Delta);
end%function

function [r, nsquare] = expmnosquare(A)
    % [r, nsquare] = expmnosquare(A)
    %
    % Returns r and nsquare such that expm(A) = r^(2^nsquare). This is useful
    % if you have recursive formulas to preform the squaring yourself, e.g.,
    % as in Computing Integrals Involving the Matrix Exponential"
    % (Van Loan, 1987).
    %
    % Base implementation is taken from Octave but without trace rescaling or
    % balancing.
    narginchk(1, 1);
    n = size(A, 1);
    
    % Scaling.
    [~, e] = log2(norm(A, 'inf'));
    nsquare = min(max(0, e), 1023);
    a = A*2^(-nsquare);

    % Pade approximation for exp(A).
    c = [5.0000000000000000e-1,...
         1.1666666666666667e-1,...
         1.6666666666666667e-2,...
         1.6025641025641026e-3,...
         1.0683760683760684e-4,...
         4.8562548562548563e-6,...
         1.3875013875013875e-7,...
         1.9270852604185938e-9];

    a2 = a^2;
    id = eye(n);
    x = (((c(8)*a2 + c(6)*id)*a2 + c(4)*id)*a2 + c(2)*id)*a2 + id;
    y = (((c(7)*a2 + c(5)*id)*a2 + c(3)*id)*a2 + c(1)*id)*a;

    r = (x - y)\(x + y);
end%function

function [allokay, args] = runtest(Delta, testtype, verbose)
    % Runs a test of the various c2d functions. testtype can be 'tricky' for
    % a particularly tricky matrix, or 'rand' for all random matrices. 'm' can
    % also be appended to the testtype to include a nonzero continuous-time
    % m matrix, although the old methods will not be tested. Finally, testtype
    % can be a struct of matrices to use directly.
    narginchk(1, 3);
    if nargin() < 2
        testtype = 'tricky';
    end
    if nargin() < 3
        verbose = false();
    end
    includem = (testtype(end) == 'm');
    if includem
        testtype = testtype(1:(end - 1));
    end
    args = struct('Delta', Delta);
    if isstruct(testtype)
        args.a = testtype.a;
        args.b = testtype.b;
        args.q = testtype.q;
        args.r = testtype.r;
        args.m = mpctools.structget(testtype, 'm', []);
    else
        switch testtype
        case 'tricky'
        args.a = [-10, 1, 0; 0, -5, 1; 0, 0, 0];
        args.b = [1, 0; 0, 0; 0, 1];
        Nx = 3;
        Nu = 2;
        case 'rand'
            Nx = randi(10);
            Nu = randi(10);
            t = rand(Nx, Nx);
            d = diag(-rand(Nx, 1));
            args.a = t\d*t;
            args.b = rand(Nx, Nu);
        otherwise
            error('Unknown testtype!');
        end
        q = rand(Nx, Nx);
        args.q = q + q';
        r = rand(Nu, Nu);
        args.r = r + r';
        if includem
            args.m = rand(Nx, Nu);
        end
    end
    
    dt = struct();
    dt.expm = mpctools.c2d('**', args);
    dt.quad = mpctools.c2d('quad', true(), 'Nquad', 1000, '**', args);
    if ~includem
        if ~isequal(testtype, 'tricky')
            dt.old = oldc2d(args);
        end
        dt.lessold = lessoldc2d(args);
    end
    allokay = checksame(dt, verbose, 1e-6);
end%function

function allokay = checksame(s, verbose, tol)
    % Checks the pairwise agreement between fields of struct s.
    narginchk(1, 3);
    if nargin() < 2
        verbose = false();
    end
    if nargin() < 3
        tol = 1e-7;
    end
    allokay = true();
    fields = fieldnames(s);
    for i = 1:length(fields)
        fi = fields{i};
        for j = (i + 1):length(fields)
            fj = fields{j};
            [ds, maxerr] = structdiff(s.(fi), s.(fj));
            okay = maxerr < tol;
            if verbose || ~okay
                fprintf('   %s / %s: %g maximum error\n', fi, fj, maxerr);
                if verbose && ~okay
                    disp(ds);
                end
            end
            allokay = allokay && okay;
        end
    end
end%function

function [ds, maxerr] = structdiff(s1, s2)
    % Computest the difference between fields of s1 and s2.
    ds = struct();
    fields = intersect(fieldnames(s1), fieldnames(s2));
    for i = 1:length(fields)
        f = fields{i};
        ds.(f) = s1.(f) - s2.(f);
    end
    if nargout() >= 2
        maxerr = max(structfun(@(x) max(abs(x(:))), ds));
    end
end%function

% Run tricky test.
fprintf('*** Tricky test ***\n');
assert(runtest(10, 'tricky', true()));

% Run random tests.
rand('state', 1000);
Ntests = 20;
okay = true(Ntests, 1);
failures = cell(Ntests, 1);
for i = 1:Ntests
    fprintf('*** Random test %d\n', i);
    if i > Ntests/4
        Delta = 50;
        testtype = 'randm';
    else
        Delta = 10;
        testtype = 'rand';
    end
    [okay(i), mats] = runtest(Delta, testtype, false());
    if ~okay(i)
        failures{i} = mats;
    end
end

