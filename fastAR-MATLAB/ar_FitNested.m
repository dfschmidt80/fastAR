%% V1.0
%% Fitted a nested sequence of AR models using the fastAR algorithm
% Based on the algorithm in the paper 
%   "Extremely Fast Maximum Likelihood Estimation of High-Order Autoregressive Models"
%   D. F. Schmidt and E. Makalic, Journal of Time Series Analysis (to appear)
% 
% function [phi, rho, psi, sigma2, mu, L] = ar_FitNested(y, pmax, varargin)
%
% Parameters:
%   y          = time series for fitting. Can be a matrix, in which columns
%                are the different time series. In this case each series
%                must be the same length.
%                Can also be a cell array, in which each series is in a one
%                of the cells; in this case, the different time series can
%                have different lengths.
%   p          = maximum order of AR model to estimate
%   varargin   = optional arguments described below.
%
% The following optional arguments are supported in the format 'argument', value:
%  'demean'    - the type of demeaning required if no 'mu' is passed: can
%                be 'off' for no mean adjustment, 'common' for a single mean
%                for all series, and 'series' for a seperate mean for each
%                series. Default: 'common'
%  'meanest'   - type of estimator for estimator 'mu' if not passed; can be
%                'sample' to use the sample mean, or 'ML' to use exact
%                maximum likelihood. Default: '
%  'reltol'    - relative tolerance of the fitting algorithm. Default: 0.005
%       
% Returns:
%   rv         = an object containing all fitted models. Can be used in
%                conjunction with ar_Select() to perform model selection
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
function rv = ar_FitNested(y, pmax, varargin)

%% Parse input
inParser = parseInput(pmax, varargin);

demean = inParser.Results.demean;
meanest = inParser.Results.meanest;
reltol = inParser.Results.reltol;

%% Pre-compute the 'D' matrices
D = ar_ProcessY(y, pmax, demean, meanest);
m = D.m;

%% Initialise
phi    = cell(1, pmax+1);
rho    = cell(1, pmax+1);
L      = zeros(1, pmax+1);
if (strcmp(demean,'series'))
    mu = zeros(m, pmax+1);
else
    mu = zeros(1, pmax+1);
end
sigma2 = zeros(1, pmax+1);

%% Fit the empty model
[~, ~, ~, sigma2(1), mu(:,1), L(1)] = ar_FitAR_core(y, pmax, demean, meanest, [], reltol, [], D);
phi{1} = [1, zeros(1, pmax)];
rho{1} = zeros(1, pmax);

%% Fit the sequence of nested models iteratively
for i = 1:pmax
    [phi{i+1}, rho{i+1}, ~, sigma2(i+1), mu(:,i+1), L(i+1)] = ar_FitAR_core(y, pmax, demean, meanest, 1:i, reltol, rho{i}, D);
end

% Store results in a return structure
rv.pmax = pmax;
rv.ntotal = sum(D.ns);
rv.type = "nested";
rv.demean = demean;
rv.phi = phi;
rv.rho = rho;
rv.mu = mu;
rv.sigma2 = sigma2;
rv.L = L;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse input 
function inParser = parseInput(pmax, args)

%% data
defaultDemean = 'common';
defaultMeanEst = 'ML';
defaultReltol = 0.005;

inParser = inputParser;    
inParser.KeepUnmatched = true;

expectedDemean = {'off', 'common', 'series'};
expectedMeanEst = {'sample', 'ML'};

addParameter(inParser,'demean', defaultDemean, @(x)any(validatestring(x,expectedDemean)));
addParameter(inParser,'meanest', defaultMeanEst, @(x)any(validatestring(x,expectedMeanEst)));
addParameter(inParser,'reltol', defaultReltol, @(x)isnumeric(x));

% parse options 
parse(inParser, args{:});    

%% done
end