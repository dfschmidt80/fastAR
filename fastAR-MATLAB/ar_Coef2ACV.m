%% V1.0
%% Generate theoretical autocovariances for a given AR(p) autoregressive model
% 
% function [gamma, d_gamma] = ar_ACV(phi, sigma2, ng)
%
% Parameters:
%   phi        = autoregressive coefficients (with leading 1) 
%   sigma2     = innovation variance 
%   ng         = number of autocovariances to generate (from gamma_0...gamma_ng) 
%
% Returns:
%   gamma      = theoretical autocovariances 
%   d_gamma    = derivatives of autocovariances wrt. to entries of 'phi' 
%                row 'i' contains the ng derivatives of gamma_0...gamma_ng
%                wrt to phi(i)
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function [gamma, d_gamma] = ar_Coef2ACV(phi, sigma2, ng)

if (phi(1) ~= 1)
    error('Leading coefficient should be a 1');
end

phi = ar_makerow(phi);
phi = phi(2:end);
p = length(phi);

%% First p elements of AR autocovariance
A1 = triu(toeplitz([1 phi]));
A2 = fliplr(tril(toeplitz(fliplr([1 phi]))));

gamma = zeros(1,ng);

gamma(1:p+1) = fliplr(((A1+A2) \ [zeros(1, p), 1]')');
gamma(1) = gamma(1)*2;

%% Next values of gamma computed recursively 
for i = (p+2):ng
    gamma(i) = - phi*fliplr(gamma(i-p:i-1))';
end

%% Compute partial derivatives w.r.t. all coefficients
if (nargout > 1)
    d_gamma = zeros(p, ng);
    
    % For each coefficient
    iA = (A1+A2)\eye(p+1);
    for k = 1:p
        %% Initial 'p+1' partial derivatives (for gamma_0, ..., gamma_p)
        j = 0:p;
        V = -gamma(abs(p - k - j) + 1);
        d_gamma(k, 1:p+1) = fliplr((iA * V')');
        d_gamma(k, 1) = d_gamma(k, 1)*2;
        
        %% Generate remaining partial derivatives recursively
        for n = (p+2):ng
            d_gamma(k, n) = -gamma(n-k);
            for i = 1:p
                d_gamma(k, n) = d_gamma(k, n) - phi(i)*d_gamma(k,n-i);
            end
        end
    end
    
    %% Scale derivatives by innovation variance
    d_gamma = d_gamma*sigma2;
    d_gamma = d_gamma(:,1:ng);
end

%% Scale autocovariances by innovation variance
gamma = gamma*sigma2;
gamma = gamma(1:ng);

return;