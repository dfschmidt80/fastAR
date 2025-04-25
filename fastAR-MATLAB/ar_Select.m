%% V1.0
%% Select an AR order using an information criterion
% 
% function [phi, sigma2, mu, L, rho, score] = ar_Select(rv, criterion)
%
% Parameters:
%   rv         = list of nested autoregressive models (produced by a 
%                call to ar_FitNested)
%   criterion  = model selection criterion to use ('AIC', 'AICc', 'BIC',
%                'KIC', 'KICc', 'CIC')
%   average    = if false, the code will return the single best model (as
%                determined by the chosen criteria)
%                if true, the code will return a model that is formed by
%                weighted averaging of all the models, weighted by
%                exp(-score) for each model. This will usually be more complex
%                (its order will be equal to the maximum considered order) but
%                will frequently predict better. Default: false.
%
% Returns:
%   phi        = estimated autoregressive model coefficients (including leading 1)
%   sigma2     = estimated innovation variance 
%   mu         = estimated mean of the time series
%   L          = negative log-likelihood of the selected model
%   rho        = estimated partial autocorrelations
%   score      = model selection criterion score of selected model
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function [phi, sigma2, mu, L, rho, score] = ar_Select(rv, criterion, average)

if (~exist('average','var'))
    average = false;
end

pmax = length(rv.rho)-1;

%% Initialise
score = array2table(zeros(pmax+1,8));
score{:,1} = (0:pmax)';
score.Properties.VariableNames = {'Order','AIC','AICc','BIC','KIC','KICc','CIC','NML'};

L = rv.L';
n = rv.ntotal;
p = (0:pmax)';

%% Fill in the scores
score{:,2} = L + p;
score{:,3} = L + n*p./(n-p-2);
score{:,4} = L + (1/2)*log(n)*p;
score{:,5} = L + p*3/2;
score{:,6} = L + (p+1)*n./(n-p-2) - n/2*psi((n-p)/2) + n/2*log(n/2);

% CIC
i = p;
i(1) = 1;
vi = 1./(n+1-i);
if (~rv.demean)
    vi(1) = 0;
end

fic_tail = 3*cumsum(vi);
fsic_tail = cumprod((1+vi)./(1-vi)) - 1;
cic_tail = max(fic_tail, fsic_tail);
score{:,7} = L + cic_tail*(n/2);

% NML
for i = 0:pmax
    if (i == 0)
        score{i+1,8} = L(i+1);
    else
        rho_max = max(abs(rv.rho{i+1}));
        score{i+1,8} = L(i+1) + p(i+1)/2*log(n/2/pi) + ceil(p(i+1)/2)*log(asin(rho_max)) + floor(p(i+1)/2)*log(atanh(rho_max)) + p(i+1)*log(2) + (1/2)*log(n);
    end
end

%% Select the model with the requested criterion
if (~average)
    [~,I] = min(score{:,criterion});

    phi = rv.phi{I}(1:I);
    mu  = rv.mu(:,I);
    sigma2 = rv.sigma2(I);
    L = rv.L(I);
    rho = rv.rho{I}(1:(I-1));

%% Weighted model averaging    
else
    s = score{:,criterion};
    s = s-min(s);
    s = exp(-s)/sum(exp(-s));
    
    % Weight the coefficients
    phi = zeros(1, pmax+1);
    for i = 1:length(s)
        phi_est = rv.phi{i};
        phi(1:length(phi_est)) = phi(1:length(phi_est)) + phi_est*s(i);
    end
    rho = [];
    
    % Weight the variances and means
    sigma = 0;
    mu = zeros(size(rv.mu,1),1);
    for i = 1:length(s)
        mu = mu + s(i)*rv.mu(:,i);
        sigma = sigma + s(i)*sqrt(rv.sigma2(i));
    end
    
    sigma2 = sigma^2;
    phi(1)=1;
end

end