%% V1.0
%% Compute theoretical autocovariances of an AR process from partial autocorrelations
% 
% function gamma = ar_PAC2ACV(rho, sigma2, ng)
%
% Parameters:
%   rho        = partial autocorrelations of AR process
%   sigma2     = innovation variance
%
% Returns:
%   ng         = number of autocovariances (starting from gamma_0)
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function gamma = ar_PAC2ACV(rho, sigma2, ng)

p = length(rho);
phi = zeros(1,p);
ns = min(p, ng-1);

gamma = zeros(1,ng);
gamma(1) = 1;

%% Compute first p+1 autocorrelations
for j = 1:ns
    % Base case
    if (j == 1)
        phi(1) = -rho(1);
    else
        a = fliplr(phi(1:j-1));
        phi(1:j-1) = phi(1:j-1) - a*rho(j);
        phi(j) = -rho(j);
    end
    
    % Compute next autocorrelation
    i = 1:j;
    d = abs(i-j)+1;
    gamma(j+1) = -gamma(d)*phi(i)';
end

%% Compute remaining autocorrelations
if (ng > (p+1))
    %% Next values of gamma computed recursively 
    for i = (p+2):ng
        gamma(i) = -phi*fliplr(gamma(i-p:i-1))';
    end
end

%% Scale by marginal variance to convert to autocovariances
g0 = exp(-sum(log(1-rho.^2)))*sigma2;
gamma = gamma*g0;
    
end