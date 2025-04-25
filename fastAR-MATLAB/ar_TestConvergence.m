%% Generate experimental results for convergence tests for Figure 2 of the paper
clear;

rng(3154);
niter = 2;

% Sample sizes to test
n   = ceil(logspace(2,log10(30000),100));
SNR = 5;

ME_sqrt_n     = zeros(niter, length(n));
ME_2_sqrt_n   = zeros(niter, length(n));
ME_n_root_2_3 = zeros(niter, length(n));

% Generate a model and control for SNR
rho = rand(1,80)*2 - 1;
rho = ar_ControlSNR(rho, sqrt(5));
phi = ar_PAC2Coef(rho);

%% Do the tests
% We use ML to estimate the means
for j = 1:niter
    y = ar_GenerateDataFromPhi(phi, 1, max(n));
    
    for i = 1:length(n)
        pmax = ceil(sqrt(n(i)));
        fprintf('%d, ', pmax);

        % Fit sqrt-n
        [phi_hat, ~, ~, ~, ~, ~] = ar_FitAR_core(y(1:n(i)), pmax, 'common', 'ml', 1:pmax, 0.005, [], []);
        ME_sqrt_n(j,i) = ar_ME(phi_hat, phi);

        % Fit 2*sqrt-n
        pmax = 2*ceil(sqrt(n(i)));
        [phi_hat, ~, ~, ~, ~, ~] = ar_FitAR_core(y(1:n(i)), pmax, 'common', 'ml', 1:pmax, 0.005, [], []);
        ME_2_sqrt_n(j,i) = ar_ME(phi_hat, phi);

        % Fit n^(2/3)
        pmax = ceil(n(i)^(2/3));
        [phi_hat, ~, ~, ~, ~, ~] = ar_FitAR_core(y(1:n(i)), pmax, 'common', 'ml', 1:pmax, 0.005, [], []);
        ME_n_root_2_3(j,i) = ar_ME(phi_hat, phi);
    end

    fprintf('\n');
end

% Save results if required
%save ar_ConvergenceTestResults.mat;