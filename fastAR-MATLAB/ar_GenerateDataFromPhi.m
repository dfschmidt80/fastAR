%% V1.0
%% Generate data from an autoregressive model
% function y = ar_GenerateDataFromPhi(phi, sigma2, n)
%
% Parameters:
%   phi        = coefficients of autoregressive model (with leading 1)
%   sigma2     = innovation variance
%   n          = length of series to be generated
%
% Returns:
%   y          = generated time series
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function y = ar_GenerateDataFromPhi(phi, sigma2, n)

%% Check
if (isempty(phi) || phi(1) ~= 1)
    error('Autoregressive coefficient vector requires the leading 1');
end

phi = ar_makerow(phi);

%% Generate data
p = length(phi);
v = normrnd(0, sqrt(sigma2), n + 1e4, 1);
y = filter(1, phi, v);
y = y(end-n+1:end);

end