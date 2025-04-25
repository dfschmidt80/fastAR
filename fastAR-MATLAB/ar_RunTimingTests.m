%% V1.0
%% Run AR timing tests (used in our paper)
% function [s_005, s_0005, s_0001, delta_L] = ar_FastARTest(p_star, p_max, n, niter)
%
% Parameters:
%   p_star     = order of true model to generate
%   p_max      = maximum order of model to try in nested search
%   n          = sample size for experiments
%   niter      = number of iterations to run the experiment for
%
% Returns:
%   s_005      = results for reltol of 0.005
%   s_0005     = results for reltol of 0.0005
%   s_0001     = results for reltol of 0.0001
%   delta_L    = average difference in negative log-likelihoods relative to the
%                reltol 0.005 model
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function [s_005, s_0005, s_0001, delta_L] = ar_RunTimingTests(p_star, p_max, n, niter, file)

y = cell(niter, 1);
if (~exist('file','var'))
    %% Create all the datasets ahead of time
    for i = 1:niter
        phi = ar_GenerateUniformCoefficients(p_star);
        y{i} = ar_GenerateDataFromPhi(phi, 1, n);
    end
else
    %% Load the data
    Y = csvread(file);
    for i = 1:niter
        y{i} = Y(i,:)';
    end
end
L = cell(1,3);
L{1} = zeros(niter, p_max+1);
L{2} = zeros(niter, p_max+1);
L{3} = zeros(niter, p_max+1);

%% Test reltol = 0.005
tic;
for i = 1:niter
    rv = ar_FitNested(y{i}, p_max, 'reltol', 0.005);
    L{1}(i,:) = rv.L;
end
s_005 = toc;

%% Test reltol = 0.0005
tic;
for i = 1:niter
    rv = ar_FitNested(y{i}, p_max, 'reltol', 0.0005);
    L{2}(i,:) = rv.L;
end
s_0005 = toc;

%% Test reltol = 0.0001
tic;
for i = 1:niter
    rv = ar_FitNested(y{i}, p_max, 'reltol', 0.0001);
    L{3}(i,:) = rv.L;
end
s_0001 = toc;

%% Delta-L
delta_L(1) = 0;
delta_L(2) = mean(mean(L{1} - L{2}));
delta_L(3) = mean(mean(L{1} - L{3}));

end