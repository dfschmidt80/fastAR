%% V1.0
%% Produce forecasts for an AR model
% 
% function yhat = ar_ForecastAR(y, phi, mu, nf, sigma2, ns)
%
% Parameters:
%   y          = time series to forecast from
%   phi        = autoregressive model coefficients (including leading 1)
%   mu         = estimated mean of the time series
%   nf         = forecast horizon
%   sigma2     = innovation variance 
%   ns         = number of realisations of the forecast to produce
%                if ns = 0 then only the mean prediction is returned
%                if ns >= 1 then the realisations are returned as well
%                (default: 0)
%
% Returns:
%   yhat_mu    = mean forecast
%   ys         = if ns >= 1, then this is an ns x nf matrix of future
%                random realisations from the AR model; these can be used
%                to find Monte-Carlo estimates of intervals,
%                transformations of the forecasts, etc.
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function [yhat_mu, ys] = ar_ForecastAR(y, phi, mu, nf, sigma2, ns)

y = ar_makerow(y);

if (~exist('ns','var'))
    ns = 0;
end

%% Generate innovations to use if variance is supplied
if (ns > 0)
    v = [zeros(1,nf);normrnd(0,1,ns,nf) .* sqrt(sigma2)];
else
    v = zeros(1,nf);
end

%% Initialise
p = size(phi,2)-1;
state = fliplr(y(end-p+1:end));
state = repmat(state,ns+1,1);

% Forecast
yhat = zeros(ns+1, nf);
for j = 1:nf
    % Forecast next value 
    yhat(:,j) = -sum(phi(2:end).*state,2) + v(:,j);
    
    % Update state
    state(:,2:end) = state(:,1:end-1);
    state(:,1) = yhat(:,j);
end
yhat = yhat + mu;

yhat_mu = yhat(1,:);
if (ns > 0)
    ys = yhat(2:end,:);
else
    ys = [];
end

yhat_mu = yhat_mu';
ys = ys';

end