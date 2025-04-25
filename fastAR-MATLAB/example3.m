%% Example 3 -- Fitting multiple models using pooling over multiple series
fprintf('Example 3====================================================================\n');

rng(10);
clear;
p_star = 20;

%% Eg. 3.1 -- equal length series
% Generate three series with (n=120) from a high order AR model
Y = zeros(120,3);
phi_true = ar_GenerateUniformCoefficients(p_star);
Y(:,1) = ar_GenerateDataFromPhi(phi_true,1,120);
Y(:,2) = ar_GenerateDataFromPhi(phi_true,1,120);
Y(:,3) = ar_GenerateDataFromPhi(phi_true,1,120);
% Store them column-wise in

fprintf('Example 3.1------------------------------------------------------------------\n');
fprintf('Multiple series (equal length) with common mean; common mean estimated via ML\n');

% Fit all models from p = 0 to p = 40 using default settings using all
% three series (i.e., pooling)
% Default setting uses a single, common mean estimated via exact ML
rv = ar_FitNested(Y, 40);

% Select using AICc
[phi, sigma2, mu, L, rho, score] = ar_Select(rv, 'AICc');
fprintf('Selected order = %d; model error = %.3g\n', length(phi)-1, ar_ME(phi_true,phi));
fprintf('\n');


%% Eg. 3.2 -- different length series
fprintf('Example 3.2------------------------------------------------------------------\n');
fprintf('Multiple series (unequal length) with common mean; common or individual means estimated via ML\n');

% If the series are different lengths we can use a cell array
clear Y;
Y{1} = ar_GenerateDataFromPhi(phi_true,1,120);
Y{2} = ar_GenerateDataFromPhi(phi_true,1,140);
Y{3} = ar_GenerateDataFromPhi(phi_true,1,80);

% Fit using defaults
rv = ar_FitNested(Y, 40);
[phi, sigma2, mu, L, rho, score] = ar_Select(rv, 'AICc');
fprintf('(Common mean via ML)      Selected order = %d; model error = %.3g\n', length(phi)-1, ar_ME(phi_true,phi));

% We could also fit using individually estimated means (one for each
% series) -- this should only be done if the series are long as the ML
% estimates of the means can be quite biased for multiple series
rv = ar_FitNested(Y, 40, 'demean', 'series');
[phi, sigma2, mu, L, rho, score] = ar_Select(rv, 'AICc');
fprintf('(Individual means via ML) Selected order = %d; model error = %.3g\n', length(phi)-1, ar_ME(phi_true,phi));
fprintf('The error gets worse as each series was generated from a model with the same mean\n\n');



%% Eg. 3.3 -- different length series with different means
fprintf('Example 3.3------------------------------------------------------------------\n');
fprintf('Multiple series (unequal length) with different means; common or individual means estimated via ML\n');

% Change the means of each series
Y{1} = Y{1} + 10;
Y{3} = Y{3} - 10;

% Fit using defaults
rv = ar_FitNested(Y, 40);
[phi, sigma2, mu, L, rho, score] = ar_Select(rv, 'AICc');
fprintf('(Common mean via ML)      Selected order = %d; model error = %.3g\n', length(phi)-1, ar_ME(phi_true,phi));

rv = ar_FitNested(Y, 40, 'demean', 'series');
[phi, sigma2, mu, L, rho, score] = ar_Select(rv, 'AICc');
fprintf('(Individual means via ML) Selected order = %d; model error = %.3g\n', length(phi)-1, ar_ME(phi_true,phi));
fprintf('The error is now better as each series was generated from a model with the same mean\n');
