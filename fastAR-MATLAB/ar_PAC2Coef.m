%% V1.0
%% Convert partial autocorrelations to autoregressive coefficients
% 
% function a = ar_PAC2Coef(rho)
%
% Parameters:
%   rho        = partial autocorrelations corresponding to the coefficients
%
% Returns:
%   a          = corresponding autoregressive coefficients (with leading 1)
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function a = ar_PAC2Coef(rho)

if (isempty(rho))
    a = 1;
    return;
end

rho = makecolumn(rho);

p = length(rho);
Y = zeros(p);
Y(1, 1) = rho(1);

% Iterate ...
for k = 2:p
    for i = 1:(k-1)
        Y(k, i) = Y(k-1, i) - rho(k)*Y(k-1, k-i);
    end
    Y(k, k) = rho(k);
end
a = [1,-Y(p,:)];

return;