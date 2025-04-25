%% V1.0
%% Fast estimation of AR models using the Maximum Likelihood
% Based on the algorithm in the paper 
%   "Extremely Fast Maximum Likelihood Estimation of High-Order Autoregressive Models"
%   D. F. Schmidt and E. Makalic, Journal of Time Series Analysis (to appear)
% 
% function [phi, rho, psi, sigma2, mu, L] = ar_FitAR(y, p, varargin)
%
% Parameters:
%   y          = time series for fitting. Can be a matrix, in which columns
%                are the different time series. In this case each series
%                must be the same length.
%                Can also be a cell array, in which each series is in a one
%                of the cells; in this case, the different time series can
%                have different lengths.
%   p          = order of AR model to estimate
%   varargin   = optional arguments described below.
%
% The following optional arguments are supported in the format 'argument', value:
%  'demean'    - the type of demeaning required if no 'mu' is passed: can
%                be 'off' for no mean adjustment, 'common' for a singl mean
%                for all series, and 'series' for a seperate mean for each
%                series. Default: 'common'
%  'meanest'   - type of estimator for estimator 'mu' if not passed; can be
%                'sample' to use the sample mean, or 'ML' to use exact
%                maximum likelihood. Default: '
%  'reltol'    - relative tolerance of the fitting algorithm. Default: 0.005
%  'rho_start' - vector of 'p' partial autocorrelations to use as a starting
%                point for the numerical search. Default: empty
%  'S'         - a vector of indices indicating which subset of partial
%                autocorrelations to estimate. The remaining partial
%                autocorrelations will be set to zero. Default is to estimate
%                all partial autocorrelations up to order 'p' (specified
%                above)
%       
% Returns:
%   L          = the negative log-likelihood of the AR model
%   mu         = estimated mean(s) of the series, if requested
%   sigma2     = estimated variance of the AR model, if requested
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function [phi, rho, psi, sigma2, mu, L] = ar_FitAR(y, p, varargin)

%% Parse input
inParser = parseInput(p, varargin);

demean = inParser.Results.demean;
meanest = inParser.Results.meanest;
reltol = inParser.Results.reltol;
rho_start = inParser.Results.start;
S = inParser.Results.subset;
D_pre = inParser.Results.D;

%% Call the core fitting algorithm
[phi, rho, psi, sigma2, mu, L] = ar_FitAR_core(y, p, demean, meanest, S, reltol, rho_start, D_pre);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse input 
function inParser = parseInput(p, args)

%% data
defaultDemean = 'common';
defaultMeanEst = 'ML';
defaultReltol = 0.005;
defaultRhoStart = [];
defaultSubset = 1:p;
defaultD = [];
defaultX = [];

inParser = inputParser;    
inParser.KeepUnmatched = true;

expectedDemean = {'off', 'common', 'series'};
expectedMeanEst = {'sample', 'ML'};

addParameter(inParser,'demean', defaultDemean, @(x)any(validatestring(x,expectedDemean)));
addParameter(inParser,'meanest', defaultMeanEst, @(x)any(validatestring(x,expectedMeanEst)));
addParameter(inParser,'reltol',defaultReltol, @(x)isnumeric(x));
addParameter(inParser,'start',defaultRhoStart, @(x)(isnumeric(x) && length(x) == p && max(abs(x)) < 1) );
addParameter(inParser,'subset',defaultSubset, @(x)(isnumeric(x) && max(x) <= p));
addParameter(inParser,'D',defaultD, @(x)(iscell(x)));
%addParameter(inParser,'X',defaultX, @(x)(isnumeric(x)));

% % parse options 
parse(inParser, args{:});    

%% done
end