%% V1.0
%% Build the Champernowe data matrix for the likelihood function
% This is an internal function and is unlikely to be called by the user 
%
% function [D, Ds, Dn] = ar_Dmtx(y, p)
%
% Parameters:
%   y          = time series to build the matrix from, as a matrix; each
%                column is a time series, each row a time point. Set the
%                length of each column using the next argument. This allows
%                for pooling of different length time series.
%   ns         = length of each of the time series (columns); extra time
%                points are ignored.
%   p          = order of the autoregressive model
%   demean     = type of demeaning: if 'off', 'common', or 'sample'
%                then the data matrices across the series are collapsed 
%                from a tensor into a single matrix; otherwise the data 
%                matrices for each series are returned in a tensor form.
%   mu_hat     = a vector of estimated means for each series; do not
%                include if means will be estimated using the fastAR code.
%
% Returns:
%   D          = data matrix/tensor of matrices
%   Ds         = matrix/tensor of first order components 
%   Dn         = number of samples used in each entry of D and Ds
%                (matrix/tensor)
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function [D, Ds, Dn] = ar_Dmtx(Y, ns, p, demean, meanest, mu_hat)

m = length(ns);

if (~exist('mu_hat','var'))
    mu_hat = zeros(1,m);
end

% If a single mean passed, expand to each series (common)
if (length(mu_hat) == 1)
    mu_hat = ones(1,m)*mu_hat;
end

D = zeros(p+1,p+1,m);
Ds = zeros(p+1,p+1,m);
Dn = zeros(p+1,p+1,m);

% Build matrices in tensors as required
for i = 1:m
    [D_i, Ds_i, Dn_i] = ar_Dmtx_y(Y(1:ns(i),i) - mu_hat(i), p);

    D(:,:,i) = D_i;
    Ds(:,:,i) = Ds_i;
    Dn(:,:,i) = Dn_i;
end

% Collapse if necessary (either not estimating mean, estimating common mean
% with ML, or estimating common/series means using sample means)
if (strcmp(demean,'off') || strcmp(demean,'common') || strcmp(meanest,'sample'))
    D = sum(D,3);
    Ds = sum(Ds,3);
    Dn = sum(Dn,3);
end

end

%% Quickly build the data matrix
function [D, Ds, Dn] = ar_Dmtx_y(y, p)

n = length(y);

% %% Form the 'D' matrix -- naive version
% D  = zeros(p+1, p+1);
% Ds = zeros(p+1, p+1);
% Dn = zeros(p+1, p+1);
% for i = 1:(p+1),
%     for j = 1:(p+1),
%         k = 0:(n-i-j+1);
%         
%         D(i,j)  = y(k+i)'*y(k+j);
%         Ds(i,j) = sum(y(k+i)+y(k+j));
%         Dn(i,j) = length(k);
%     end
% end

%% Alternative algorithm
%D
D  = zeros(p+1, p+1);
D(:,1) = ar_ObservedACV(y,p)*n;
D(1,:) = D(:,1);
for j = 1:p
    for i = 1:j
        % Remove observations
        D(i+1, j+1) = D(i,j) - y(i)*y(j) - y(end-i+1)*y(end-j+1);
        D(j+1, i+1) = D(i+1, j+1);
    end
end

Ds = zeros(p+1,p+1);
ys = zeros(p+1,1);
ys(1,1) = sum(y);
for j = 2:(p+1)
    ys(j,1) = ys(j-1,1) - y(j-1) - y(end-(j-2));
end
for j = 1:(p+1)
    for i = 1:j
        % Remove observations
        Ds(i, j) = ys(i) + ys(j);
        Ds(j, i) = Ds(i, j);
    end
end

Dn = zeros(p+1,p+1);
for j = 1:(p+1)
    for i = 1:j
        Dn(i,j) = n-i-j+2;
        Dn(j,i) = Dn(i,j);
    end
end
    
end