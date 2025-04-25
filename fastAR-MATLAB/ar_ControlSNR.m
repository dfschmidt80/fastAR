%% V1.0
%% Re-scale an AR model to attain a specified signal-to-noise ratio
% function rho = ar_ControlSNR(rho, SNR)
%
% Parameters:
%   rho        = partial autocorrelations of autoregressive model
%   SNR        = target signal-to-noise ratio sqrt(E[Y^2] - 1)
%
% Returns:
%   rho        = rescaled partial autocorrelations such that
%                sqrt(prod(1./(1-rho.^2)) - 1) = SNR
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function rho = ar_ControlSNR(rho, SNR)

U = 1./max(abs(rho));
c = fminbnd(@(c)((sqrt(prod(1./(1-(c*rho).^2)) - 1) - SNR)^2), 1e-4, U - 1e-4);
rho = rho*c;

end