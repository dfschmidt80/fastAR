%% V1.0
%% Calculate partial derivatives of an autoregressive coefficient vector wrt to a partial autocorrelation
% This is the MATLAB version -- for substantial speedups please obtain the
% MEX file (or C source, and compile to MEX).
% 
% function alpha = ar_PAC_d(rho, k)
%
% Parameters:
%   rho        = partial autocorrelations corresponding to the coefficients
%   k          = which partial autocorrelation the derivatives are wrt
%
% Returns:
%   alpha      = derivatives of coefficients (sans the leading 1 term) wrt to partial autocorrelation rho_k
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function alpha = ar_PAC_d(rho, k)

p = length(rho);
alpha = zeros(1, p);

%% Initialise
if (k > 1)
    alpha(1) = -rho(1);
else
    alpha(1) = -1;
end

%% Recur ...
for j = 2:p
    if (k > j)
        alpha(1:j) = [alpha(1:j-1) - fliplr(alpha(1:j-1))*rho(j), -rho(j)];
    elseif (k == j)
        alpha(1:j) = [-fliplr(alpha(1:j-1)), -1];
    else
        alpha(1:j) = [alpha(1:j-1) - fliplr(alpha(1:j-1))*rho(j), 0];
    end
end

return;