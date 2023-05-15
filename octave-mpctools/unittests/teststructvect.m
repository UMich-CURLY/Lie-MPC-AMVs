% Tests VarLayout behavor with numeric and symbolic matrices. Makes sure it
% matches the reference implementation in structvect and istructvect.
clear('classes');
status = 0;

% Need some helper functions.
function [sym, num] = getvarsym(varargin)
    for i = 2:2:length(varargin)
        varargin{i} = {varargin{i}};
    end
    args = struct(varargin{:});
    varnames = fieldnames(args);
    n = 0;
    sym = struct();
    num = struct();
    for i = 1:length(varnames)
        % Get sizes.
        v = varnames{i};
        [vsize, vrepeat] = args.(v){:};
        
        % Get numeric values.
        vshape = [vsize, vrepeat];
        nv = prod(vshape);
        num.(v) = reshape((n + 1):(n + nv), vshape);
        n = n + nv;
        
        % Make symbols.
        if isscalar(vrepeat)
            vrepeat = [1, vrepeat];
        end
        vfun = @(t) casadi.SX.sym(sprintf('%s_%d', v, t), vsize);
        sym.(v) = arrayfun(vfun, reshape(1:prod(vrepeat), vrepeat), ...
                           'UniformOutput', false());
    end%function
end%function

function toc(name)
    % Calls toc and prints message.
    endtime = builtin('toc');
    fprintf('%20s took %.5g s.\n', name, endtime);
end%function

Nt = 100;
Nx = 10;
Nc = 5;
Nu = 15;

% All vector sequences to compare to old method.
[sym, vars] = getvarsym('x', {Nx, Nt}, 'xc', {Nx, [Nc, Nt]}, 'u', {Nu, Nt});
layout = mpctools.VarLayout(sym);
sizes = structfun(@size, vars, 'UniformOutput', false());

fprintf('\n*** Testing vec\n');

tic();
v = layout.vec(vars);
toc('VarLayout.vec')

tic();
vcheck = structvect(vars);
toc('structvect');

assert(v, vcheck);
if isequal(v, vcheck)
    disp('vec are equal.');
else
    disp('vec are NOT equal!');
    status = status + 1;
end

fprintf('\n*** Testing ivec\n')

tic();
s = layout.ivec(v);
toc('VarLayout.ivec');

tic();
scheck = istructvect(v, sizes);
toc('istructvect');

assert(s, scheck)
if isequal(s, scheck)
    disp('ivec are equal.');
else
    disp('ivec are NOT equal!');
    status = status + 1;
end

% Try with symbolics.
Nt = 25;
Nx = 10;
Na = 5; % For matrices.

[syms, nums] = getvarsym('x', {Nx, Nt}, 'A', {[Na, Nx], Nt});

layout = mpctools.VarLayout(syms);
snums = layout.vec(nums);

fprintf('\n*** Testing sym vec\n');

tic();
s = layout.symvec(syms);
toc('sym VarLayout.vec');

tic();
scheck = structvect(syms);
toc('sym structvect');

% Need to substitute numbers back in to see if things make sense.
ss = s;
sscheck = scheck;
for t = 1:Nt
    ss = casadi.substitute(ss, syms.x{t}, nums.x(:,t));
    sscheck = casadi.substitute(sscheck, syms.x{t}, nums.x(:,t));
    
    ss = casadi.substitute(ss, syms.A{t}, nums.A(:,:,t));
    sscheck = casadi.substitute(sscheck, syms.A{t}, nums.A(:,:,t));
end
ss = full(casadi.DM(ss));
sscheck = full(casadi.DM(sscheck));

assert(ss, sscheck);
if isequal(ss, sscheck)
    disp('syms are equal.');
else
    disp('syms are NOT equal!');
    status = status + 1;
end

assert(ss, snums);
if isequal(ss, snums)
    disp('sym and num are the same.');
else
    disp('sym and num are not the same.');
    status = status + 1;
end

% Check final status.
if status == 0
    fprintf('\n*** All tests successful.\n');
else
    fprintf('\n*** %d tests FAILED!\n', status);
end

