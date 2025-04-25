%% V1.0
%% Process multiple time series in preparation for model fitting
% This is an internal function and is unlikely to be called by the user 
%
% function D = ar_ProcessY(y, pmax, demean, meanest)
%
% Parameters:
%   y          = time series for fitting. Can be a matrix, in which columns
%                are the different time series. In this case each series
%                must be the same length.
%                Can also be a cell array, in which each series is in a one
%                of the cells; in this case, the different time series can
%                have different lengths.
%   pmax       = maximum order of autoregressive model that you plan to try
%   demean     = type of demeaning: 
%                 'off' for none; 
%                 'common' for a single mean for all series; 
%                 'series' for a seperate mean for each series
%   meanest    = type of estimator for the means of the series (if requested)
%                 'sample' uses the sample mean
%                 'ml' uses exact maximum likelihood estimation
%
% Returns:
%   D          = structure containing data matrices and statistics for
%                use in the fast fitting algorithm
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function D = ar_ProcessY(y, pmax, demean, meanest)

%% First process the y -- could be multiple series
% Cell array -- series of different lengths
if (iscell(y))
    m = length(y);
    
    ybar_common = 0;
    ybar_series = zeros(m,1);
    
    % First find series lengths and means
    ns = zeros(1, length(y));
    for i = 1:m
        ns(i) = length(y{i});
        ybar_series(i) = mean(y{i});
        ybar_common = ybar_common*ns(i);
    end
    
    ybar_common = ybar_common / sum(ns);
    
    % Convert series into a matrix
    Y = zeros(max(ns), m);
    for i = 1:m
        Y(1:ns(i),i) = y{i};
    end
    
% Matrix -- every series same length
else
    ybar_common = mean(y(:));
    ybar_series = mean(y,1)';
    
    [n, m] = size(y);
    ns = ones(1,m)*n;
    Y = y;
end

%% Pre-compute the 'D' matrices
if (~strcmp(demean,'off') && strcmp(meanest,'sample'))
    if (strcmp(demean,'common'))
        [Dm, Ds, Dn] = ar_Dmtx(Y, ns, pmax, demean, meanest, ybar_common);
    else
        [Dm, Ds, Dn] = ar_Dmtx(Y, ns, pmax, demean, meanest, ybar_series);
    end
else
    [Dm, Ds, Dn] = ar_Dmtx(Y, ns, pmax, demean, 0);
end

%% Store
D = struct();

D.Y  = y;
D.ns = ns;
D.m  = length(ns);
D.n  = min(ns);
D.ybar_common = ybar_common;
D.ybar_series = ybar_series;

D.p = pmax;
D.D = Dm;
D.Ds = Ds;
D.Dn = Dn;

D.demean = demean;
D.meanest = meanest;

end