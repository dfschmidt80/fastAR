%% V1.0
%% Transform autoregressive coefficients into equivalent partial autocorrelations
% function rho = ar_CoefficientsToPAC(phi)
%
% Parameters:
%   phi        = autoregressive coefficients (with leading 1)
%
% Returns:
%   rho        = partial autocorrelations of autoregressive coefficients
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function rho = ar_Coef2PAC(phi)

if (phi(1) ~= 1)
    error('Leading coefficient should be a 1');
end

phi = ar_makerow(phi);
phi = phi(2:end);
p = length(phi);

% Compute the first (p+1) autocorrelations
gamma = ar_Coef2ACV([1,phi], 1, p+1);
autocorr = gamma / gamma(1);

N = zeros(p);
rho = zeros(1,p);

% Now compute the PACFs
for i = 1:p
    % Build the numerator matrix
    for j = 1:i
        lag = 0 - (j-1);
        for k = 1:i
            N(k, j) = autocorr(abs(lag) + 1);
            lag = lag + 1;
        end
    end
    
    % Solve for PACF
    M = N(1:i,1:i)\autocorr(2:i+1)';
    rho(i) = M(end);
end

return;