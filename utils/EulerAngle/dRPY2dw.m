function wb = dRPY2dw(rpy, drpy)
x = rpy(1);
y = rpy(2);
z = rpy(3);

wx = drpy(1);
wy = drpy(2);
wz = drpy(3);

wb = [wx - wz*sin(y);...
      wy*cos(x) + wz*cos(y)*sin(x);...
      wz*cos(x)*cos(y) - wy*sin(x)];
end