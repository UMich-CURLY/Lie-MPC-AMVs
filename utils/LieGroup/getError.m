function [eX, ex, exi] = getError(X, X0, xi, xi0)
%% X^{-1} * X0
eX = X^(-1) * X0;
ex = logm(eX);
exi = xi - xi0;
end