function adx = adjoint(x)
w = [x(1);x(2);x(3)];
v = [x(4);x(5);x(6)];
adx = [skew(w), zeros(3);...
       skew(v), skew(w)];
end
