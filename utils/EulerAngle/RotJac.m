%%
syms x y z 'real'
syms wx wy wz 'real'

R = rotz(z)  * roty(y) * rotx(x);
dR = R;
for i = 1:3
    for j = 1:3
        dR(i, j) = jacobian(R(i, j), [x, y, z]) * [wx, wy, wz]';
    end
end
ww = R' * dR;
ww = simplify(ww)
ww(3,2)
ww(1,3)
ww(2,1)
%%
function R = rotx(t)
R = [1, 0, 0;...
     0, cos(t), -sin(t);...
     0, sin(t), cos(t)];
end

function R = roty(t)
R = [cos(t), 0, sin(t);...
     0, 1, 0;...
     -sin(t), 0, cos(t)];
end

function R = rotz(t)
R = [cos(t), -sin(t), 0;...
     sin(t), cos(t), 0;...
     0,0,1];
end