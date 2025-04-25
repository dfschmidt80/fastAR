%% Example 1 -- Fitting a single model
%% Eg. 1.1 -- different tolerances and nested fitting
fprintf('Example 1====================================================================\n');
fprintf('Ex. 1.1----------------------------------------------------------------------\n');
rng(1);
clear;
p_star = 40;

fprintf('** True model coefficients\n');
% Generate an n=500 time series from a high order AR model
phi_true = ar_GenerateUniformCoefficients(p_star)
y = ar_GenerateDataFromPhi(phi_true,1,5e2);
pause

% Fit just a p_star-th order model, varying whether a mean is estimated, and
% the relative tolerance
fprintf('** Estimated coefficients at rel_tol = 0.005 with mean estiamted via ML\n');
[phi, ~, ~, ~, mu, L] = ar_FitAR(y,p_star)
pause

fprintf('** Estimated coefficients at rel_tol = 0.005 with no mean (i.e., mu = 0)\n');
[phi, ~, ~, ~, mu, L] = ar_FitAR(y,p_star,'demean','off')
pause

fprintf('** Estimated coefficients at rel_tol = 0.0001 with mean estiamted via ML\n');
[phi, ~, ~, ~, mu, L] = ar_FitAR(y,p_star,'demean','off','reltol',1e-4)
fprintf('Note: decreasing reltol can improve the negative log-likelihood,\n');
fprintf('particularly if you are fitting a single model with a cold-start\n\n');
pause

fprintf('** We could also fit a p=40 model as part of a nested sequence; here we fit all\n');
fprintf('** models from p = 0 to p = 60, and look at the p=40 model\n');

% Fit all models from 0 through to p_star+20 with no mean using default relative
% tolerance
tic;rv = ar_FitNested(y, p_star+20, 'demean', 'off');toc
% Takes 0.05 seconds on my machine to fit 61 AR models ...

% Let's examine the negative log-likelihood of the AR(40) model estimated
% in the sequence
rv.L(p_star+1)
fprintf('Note the negative log-likelihood of the p=40 model in the sequence, using default reltol is\n'); 
fprintf('better than fitting a single model at the same reltol (and almost as good as the much smaller\n');
fprintf('reltol used above) -- this can sometimes happen due to the warm starting in the\n'); 
fprintf('nested sequence leading to better starting points\n\n');

pause

%% Least-squares fit for comparison
%% Eg. 1.2 -- ML vs least-squares
fprintf('Ex. 1.2----------------------------------------------------------------------\n');

fprintf('** Fitting the model via least-squares\n');
Q=lagmatrix(y,0:p_star);
phi_ls = (Q(p_star+1:end,2:end)\Q(p_star+1:end,1))'

fprintf('** Comparison of ML and LS\n');
fprintf('Model error for ML = %.3g\n', ar_ME(phi_true,phi));
fprintf('Model error for LS = %.3g\n', ar_ME(phi_true,[1,-phi_ls]));

fprintf('\nIn this case ML has a %.3g%% smaller error than least-squares\n', 100*(1-(ar_ME(phi_true,phi) / ar_ME(phi_true,[1,-phi_ls]))));