function J = Jac_SO3_L(xi)
t = norm(xi);
J = eye(3) + (1 - cos(t)) / t^2 * skew(xi) + (t - sin(t)) / t^3 * skew(xi)^2;
end