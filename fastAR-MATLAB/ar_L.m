%% V1.0
%% Compute the negative log-likelihood for an AR model
% 
% function [L, mu, sigma2] = ar_L(y, theta, varargin)
%
% Parameters:
%   y          = time series for fitting. Can be a matrix, in which columns
%                are the different time series. In this case each series
%                must be the same length.
%                Can also be a cell array, in which each series is in a one
%                of the cells; in this case, the different time series can
%                have different lengths.
%   theta      = autoregressive model parameters; either coefficients with
%                leading 1 (default), or partial autocorrelations
%   varargin   = optional arguments described below.
%
% The following optional arguments are supported in the format 'argument', value:
%  'thetatype' - how to interpret 'theta'; this can be: 'phi' if theta
%                are autoregressive coefficients (including the lead 1), or
%                'rho' if theta are the partial autocorrelations
%  'mu'        - series mean(s). Must either be a constant if a common mean
%                is used for all series, or a vector with as many entries
%                as the number of series. If it is not passed, the mean(s)
%                will be estimated from the data 'y'
%  'sigma2'    - the innovation variance. If this is not passed, it will be
%                estimated from the data.
%  'demean'    - the type of demeaning required if no 'mu' is passed: can
%                be 'off' for no mean adjustment, 'common' for a singl mean
%                for all series, and 'series' for a seperate mean for each
%                series. Default: 'common'
%  'meanest'   - type of estimator for estimator 'mu' if not passed; can be
%                'sample' to use the sample mean, or 'ML' to use exact
%                maximum likelihood. Default: 'ML'
%
% Returns:
%   L          = the negative log-likelihood of the AR model
%   mu         = estimated mean(s) of the series, if requested
%   sigma2     = estimated variance of the AR model, if requested
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function [L, mu, sigma2] = ar_L(y, theta, varargin)

%% Parse input
inParser = parseInput(varargin);

theta_type = inParser.Results.thetatype;
demean = inParser.Results.demean;
meanest = inParser.Results.meanest;
sigma2 = inParser.Results.sigma2;
mu = inParser.Results.mu;

%% If partial autocorrelations are passed
if (strcmp(theta_type,'rho'))
    rho = theta;
    phi = ar_PAC2Coef(rho)';
else
    rho = ar_Coef2PAC(theta);
    phi = ar_makecolumn(theta);
end

% Prepare data
p = length(phi)-1;

D = ar_ProcessY(y, p, demean, 'ML');
m = D.m;
n = sum(D.ns);

if (~isempty(mu) && length(mu) ~= 1 && length(mu) ~= m)
    error('If a mean vector is passed it must either be a scalar (common mean for all series) or a vector (one for each series).');
end

% Estimate 'mu', if not passed
[mu, Dstar] = makeDstar(phi, D, mu, demean, meanest);

% Estimate 'sigma2' if not passed
S1 = phi'*Dstar*phi;
if (isempty(sigma2))
    sigma2 = S1/n;
end

% Likelihood
L = (n/2)*log(2*pi*sigma2) - (m/2)*sum( (1:p) .* log(1-rho.^2) ) + S1/sigma2/2;

end

%% Build the D-star matrix based on mean estimation options
function [mu_hat, Dstar] = makeDstar(beta, D, mu, demean, meanest)

% Use sample means if requested
if (strcmp(demean,'off') || (isempty(mu) && strcmp(meanest,'sample')))
    if (strcmp(demean,'off'))
        mu = 0;
        demean = 'common';
    else
        if (strcmp(demean,'common'))
            mu = D.ybar_common;
        else
            mu = D.ybar_series;
        end
    end
end

%% Common
if (strcmp(demean,'common'))
    % Estimate using ML, if required
    if (isempty(mu))
        mu_hat = beta'*D.Ds*beta/2/(beta'*D.Dn*beta);
    else
        mu_hat = mu;
    end
    Dstar = D.D - mu_hat*D.Ds + mu_hat^2*D.Dn;

%% Series
elseif (strcmp(demean,'series'))
    % Estimate using ML, if required
    mu_hat = zeros(1, 1, D.m);
    if (isempty(mu))
        mu_hat = zeros(1, 1, D.m);

        num = pagemtimes(pagemtimes(beta',D.Ds),beta);
        den = pagemtimes(pagemtimes(beta',D.Dn),beta);
        mu_hat(1,1,:) = (num(:)./den(:)/2)';
    else
        mu_hat(1,1,:) = mu;
    end

    Dstar = D.D - pagemtimes(mu_hat,D.Ds) + pagemtimes((mu_hat.^2),D.Dn);
    Dstar = sum(Dstar,3);

    % Convert to vector
    mu_hat = mu_hat(:);
    
%% No demeaning
else
    Dstar = D.D;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse input 
function inParser = parseInput(args)

%% data
defaultThetatype = 'phi';
defaultMu = [];
defaultDemean = 'common';
defaultMeanEst = 'ML';
defaultSigma2 = [];
%defaultD = [];

inParser = inputParser;    
inParser.KeepUnmatched = true;

expectedThetatype = {'phi', 'rho'};
expectedDemean = {'off', 'common', 'series'};
expectedMeanEst = {'sample', 'ML'};

addParameter(inParser,'thetatype', defaultThetatype, @(x)any(validatestring(x,expectedThetatype)));
addParameter(inParser,'demean', defaultDemean, @(x)any(validatestring(x,expectedDemean)));
addParameter(inParser,'meanest', defaultMeanEst, @(x)any(validatestring(x,expectedMeanEst)));
addParameter(inParser,'sigma2', defaultSigma2, @(x)(isnumeric(x)));
addParameter(inParser,'mu', defaultMu, @(x)(isnumeric(x)));

% % parse options 
parse(inParser, args{:});    

%% done
end