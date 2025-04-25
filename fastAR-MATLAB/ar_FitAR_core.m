%% V1.0
%% Fast estimation of AR models using the Maximum Likelihood core fitting code
% This is an internal function with limited error checking -- users should generally use
% the ar_FastAR() function instead to fit a single model.
% 
% function [phi, rho, psi, sigma2, mu, L] = ar_FitAR_core(y, p, demean_str, meanest_str, subset, reltol, start, D_pre)
%
% Parameters:
%   y           = time series to estimate model from
%   p           = order of model to estimate
%   demean_str  = {'off', 'common', 'series'}
%   meanest_str = {'sample', 'ML'}
%   subset      = vector with which partial autocorrelations to estimate
%   reltol      = relative tolerance for stopping fitting
%   start       = vector of partial autocorrelations to start from
%   D_pre       = pre-computed D matrices (if [], will be generated from
%                 the data)
%
% Returns:
%   rho         = estimated partial autocorrelations at the posterior mode
%   tau         = estimated innovation variance at the posterior mode
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function [phi, rho, psi, sigma2, mu, L] = ar_FitAR_core(y, p, demean_str, meanest_str, subset, reltol, start, D_pre)

%% Process the input time series
%[Y, ns, ybar_common, ybar_series] = ar_ProcessY(y);
%n_min = min(ns);
%ntotal = sum(ns);
%m = length(ns);

%% Construct data matrix, if vector of data passed
if (isempty(D_pre))
    D = ar_ProcessY(y, p, demean_str, meanest_str);
else
    %% Else we use the provided data
    D = D_pre;
end

n_min = min(D.ns);
ntotal = sum(D.ns);
m = D.m;

%% Argument checking
pmax = floor(n_min/2);
if (~isnumeric(p) || p < 0 || p > pmax)
    error('Order p must be between 0 and floor(n/2)');
end

rho_start = start;
S = subset;

%% If a single lambda has been passed
psi = [];

%% Begin by fitting the ML model, if no starting point passed
if (isempty(rho_start))
    rho = zeros(1,p);
else
    %% Else use provided starting point
    rho = rho_start;
end

% %% Begin by fitting regression terms using least-squares if no starting point passed
% if (~isempty(X))
%     if (~exist('psi_start','var'))
%         psi = X\(y-mean(y));
%     else
%         %% Else use provided starting point
%         psi = psi_start;
%     end
%     %[D, Ds, Dn] = ar_Dmtx(y - X*psi, p);
% end
    
%% Process demean/meanest
demean = 0;
if (strcmp(demean_str,'common'))
    demean = 1;
elseif (strcmp(demean_str,'series'))
    demean = 2;
end
meanest = 0;
if (strcmp(meanest_str,'sample'))
    meanest = 1;
end

%% Initial estimate for the mean and variance
beta  = ar_PAC2Coef(rho)';
[mu, Dstar] = estimateMean(beta, D, demean, meanest);
tau   = beta'*Dstar*beta/ntotal;

%% If empty model is requested, we are done
if (p == 0)
    % Means -- common or series    
    %if (demean == 1)
    %    mu = ybar_common;
    %else
    %    mu = ybar_series;
    %end    
    
    % Empty model stats
    rho = [];
    sigma2 = tau;
    phi = 1;
    L = ntotal/2*log(2*pi*sigma2) + ntotal/2;
    
    return;
end

%% Co-ordinate wise descent
ActiveSet     = S; 
gg            = zeros(p,1);
UseLinesearch = true;

iter = 0;
while (1)
    %% Update the rho's
    rold = rho;
    for i = ActiveSet
        % Update the beta's, if either:
        %  (i) rho(i) is not already zero, or
        %  (ii) if the previous rho(i-1) changed
        rho0 = rho;
        rho0(i) = 0;
        beta = ar_PAC2Coef(rho0)';
        
        % Calculate the derivatives of beta w.r.t. rho(i)
        g = [0;ar_PAC_d(rho0,i)'];
        
        % Compute gradient and (adjusted by L2 shrinkage) Hessian for this parameter
        gg(i) = beta'*Dstar*g/tau;
        HH    = g'*Dstar*g/tau;

        % Solve the 1-d optimisation problem by solving the cubic
        %   (HH)*r^3 + (gg)*r^2 - (HH+i)*r - (gg) = 0
        % for the single root in (-1,1)
        %
        % Only solve for the cubic if the solution will be non-zero
        C = [HH, (gg(i)), -(HH+i*m), -(gg(i))];
        [r(1),r(2),r(3)]  = cubicroots(C(2)/C(1), C(3)/C(1), C(4)/C(1));
        
        % Find the one valid root
        if (r(1) > -1 && r(1) < 1)
            rho(i) = r(1);
        elseif (r(2) > -1 && r(2) < 1)
            rho(i) = r(2);
        else
            rho(i) = r(3);
        end
    end
    
    %% Update the mean and variance
    beta = ar_PAC2Coef(rho)';
    [~, Dstar] = estimateMean(beta, D, demean, meanest);
    
    S1 = beta'*Dstar*beta;
    tau = S1/ntotal;
    
    %% Add an optional line-search
    if (1)
        rho_delta = rho - rold;
        step = [1.1,2,4,8];
        
        rho_best = rho;
        L_old = (ntotal/2)*log(2*pi*tau) - sum((1:p).*log(1-rho.^2))/2  + S1/2/tau;
        L_best = L_old;
        
        for i = 1:length(step)
            rho_new = rold + rho_delta*step(i);
            
            beta_new = ar_PAC2Coef(rho_new)';
            S1 = beta_new'*Dstar*beta_new;
            
            L_new = (ntotal/2)*log(2*pi*tau) - sum((1:p).*log(1-rho_new.^2))/2  + S1/2/tau;
            
            if (L_new < L_best)
                beta_best = beta_new;
                rho_best = rho_new;
                L_best = L_new;
            else
                break;
            end
        end        
        
        % If we updated, re-update 
        if (L_best < L_old)
            rho = rho_best;
            beta = beta_best;
            [~, Dstar] = estimateMean(beta, D, demean, meanest);

            S1 = beta'*Dstar*beta;
            tau = S1/ntotal;
        end
    end    
    
%     %% Every 2 iterations we update the regression coefficients
%     if (~isempty(X) && mod(loopcnt,5) == 0)
%         Z = filter(beta,1,[X,y-mu]);
% 
%         % Exact fit?
%         Gp = toeplitz(ar_ACVFromPAC(rho, 1, p));
%         V = Z(p+1:end,1:end-1);
%         GL = chol(Gp);
%         Q = (GL\eye(p));
%         Q = X(1:p,:)'*(Q*Q');
%         psi = (Q*X(1:p,:) + V'*V) \ (Q*(y(1:p)-mu) + V'*(Z(p+1:end,end)));
%         
%         % Update D matrices
%         [D, Ds, Dn] = ar_Dmtx(y - X*psi, p);
%     end
    
    iter = iter+1;
    %% If no significant improvement has been made, stop    
    if (sum(abs(rold-rho)) / (1+sum(abs(rho))) < reltol)
        break;
    end
end
iter;

%% We are done -- so compute final estimates
phi   = ar_PAC2Coef(rho);
beta  = phi';
[mu, Dstar] = estimateMean(beta, D, demean, meanest);

S1    = beta'*Dstar*beta;
tau   = S1/ntotal;

L     = (ntotal/2)*log(2*pi*tau) - (m/2)*sum( (1:p).*log(1-rho.^2) ) + S1/tau/2;
sigma2 = tau;

end

function [mu_hat, Dstar] = estimateMean(beta, D, demean, meanest)

% We only need to compute mean estimates if we are using ML
if (meanest == 0 && demean > 0)
    % Common
    if (demean == 1)
        mu_hat = beta'*D.Ds*beta/2/(beta'*D.Dn*beta);
        Dstar = D.D - mu_hat*D.Ds + mu_hat^2*D.Dn;
    
    % Series
    else
        mu_hat = zeros(1, 1, D.m);
        
        num = pagemtimes(pagemtimes(beta',D.Ds),beta);
        den = pagemtimes(pagemtimes(beta',D.Dn),beta);
        mu_hat(1,1,:) = (num(:)./den(:)/2)';
        
        Dstar = D.D - pagemtimes(mu_hat,D.Ds) + pagemtimes((mu_hat.^2),D.Dn);
        Dstar = sum(Dstar,3);
        
        % Convert to vector
        mu_hat = mu_hat(:);
    end
    
% Otherwise we are either not estimating the mean, or have estimated it via sample means    
else
    Dstar = D.D;
    if (demean == 0)
        mu_hat = 0;
    else
        if (demean == 1)
            mu_hat = D.ybar_common;
        else
            mu_hat = D.ybar_series;
        end
    end
end

end


%% Solve a cubic analytically
function [x1,x2,x3] = cubicroots(a, b, c)

% Determine if it is the three real roots case or one real, two complex
Q  = (a^2 - 3*b)/9;
R  = (2*a^3 - 9*a*b + 27*c)/54; 

Q3 = Q^3;
R2 = R^2;

% If three real roots, use trigonometric solutions
if (R2 < Q3)
    theta = acos(R/sqrt(Q3));
    SQ    = sqrt(Q);
    x1  = -2*SQ*cos(theta/3) - a/3;
    x2  = -2*SQ*cos((theta + 2*pi)/3) - a/3;
    x3  = -2*SQ*cos((theta - 2*pi)/3) - a/3;
    
    return
end

% Otherwise we have one real, two complex
A = -(abs(R) + sqrt(R2-Q3))^(1/3);
if (R < 0)
    A = -A; 
end
if (A==0)
    B = 0;
else
    B = Q/A;
end

AB   = A + B;
x1   = AB - a/3;
comp = 1i*sqrt(3)*(A-B)/2;
x2   = -0.5*AB - a/3 + comp;
x3   = -0.5*AB - a/3 - comp;

end