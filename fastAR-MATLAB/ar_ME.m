%% V1.0
%% Compute the (unnormalized) model error between a true autoregressive model and an estimated autoregressive model
% 
% function me = ar_ME(phi_true, phi_est, normalize)
%
% Parameters:
%   phi_true   = autoregressive coefficients of the true, generating model (with leading 1)
%   phi_est    = autoregressive coefficients of the estimated model (with leading 1)
%   normalize  = normalize the model error by dividing by marginal variance
%                of true model (default: true)
%
% Returns:
%   me         = the (unnormalized) model error
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function me = ar_ME(phi_true, phi_est, normalize)

if (phi_est(1) ~= 1 || phi_true(1) ~= 1)
    error('Leading coefficient should be a 1');
end

if (~exist('normalize','var'))
    normalize = true;
end

phi_est = ar_makerow(phi_est);
phi_true = ar_makerow(phi_true);

[phi_est,phi_true,m] = zeropadding(phi_est,phi_true);

%% Compute the model error
G = toeplitz(ar_Coef2ACV(phi_true, 1, m));
me = (phi_est - phi_true)*G*(phi_est - phi_true)';
if (normalize)
    me = me/G(1,1);
end

end

%% Pad vectors with zero to match lengths
function [an1,an2,m] = zeropadding(a1,a2)

m = max(length(a1),length(a2));
an1 = zeros(1,m);
an1(1:length(a1)) = a1;
an2 = zeros(1,m);
an2(1:length(a2)) = a2;

end