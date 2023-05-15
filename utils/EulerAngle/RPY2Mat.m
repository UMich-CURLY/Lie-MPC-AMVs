function R = RPY2Mat(rpy)
% R = rotx(rpy(1)) * roty(rpy(2)) * rotz(rpy(3));
R = rotz(rpy(3)) * roty(rpy(2)) * rotx(rpy(1));
end
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