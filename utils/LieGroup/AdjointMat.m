function Adx = AdjointMat(x)
R = x(1:3,1:3);
p = x(1:3,4);
Adx = [R, zeros(3);...
       skew(p) * R, R];
end