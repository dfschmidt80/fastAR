%% Example 2 -- Fitting multiple models and order selection
%% Eg 2.1 -- selecting models with different information criteria
fprintf('Example 2====================================================================\n');
fprintf('Eg. 2.1 ---------------------------------------------------------------------\n');

rng(1);
clear;
p_star = 20;

% Generate an n=500 time series from a high order AR model
phi_true = ar_GenerateUniformCoefficients(p_star)
y = ar_GenerateDataFromPhi(phi_true,1,5e2);

% Fit all models from p = 0 to p = 40 using default settings
% True model has p = 20
rv = ar_FitNested(y, 40);

% Select using AICc
[phi, sigma2, mu, L, rho, score] = ar_Select(rv, 'AICc');
fprintf('Selected order (AICc) = %d; model error = %.3g\n', length(phi)-1, ar_ME(phi_true,phi));

% Select using BIC
[phi, sigma2, mu, L, rho, score] = ar_Select(rv, 'AICc');
fprintf('Selected order (BIC)  = %d; model error = %.3g\n', length(phi)-1, ar_ME(phi_true,phi));

% Select using CIC
[phi, sigma2, mu, L, rho, score] = ar_Select(rv, 'CIC');
fprintf('Selected order (CIC)  = %d; model error = %.3g\n', length(phi)-1, ar_ME(phi_true,phi));

% Estimate using CIC and averaging
[phi, sigma2, mu, L, rho, score] = ar_Select(rv, 'CIC', true);
fprintf('Model averaging via CIC;    model error = %.3g\n', ar_ME(phi_true,phi));

fprintf('\n');

%% Eg 2.2 -- forecasting
fprintf('Eg. 2.2 ---------------------------------------------------------------------\n');

% First, forecast only the conditional mean
[yhat_mu, ys] = ar_ForecastAR(y, phi, mu, 10, sigma2);
fprintf('Conditional mean forecast for next 10 observations:\n');
yhat_mu'

% Get intervals via simulation
fprintf('95%% forecast interval for next 10 observations:\n');
[yhat_mu, ys] = ar_ForecastAR(y, phi, mu, 10, sigma2, 1e5);
prctile(ys', [2.5,97.5])

