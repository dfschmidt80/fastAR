%% V1.0
%% Randomly generate an autoregressive model by sampling uniformly from the stationarity region in coefficient space
% 
% function [phi, rho] = ar_GenerateARUniformCoefficients(p)
%
% Parameters:
%   p          = order of autoregressive model
%
% Returns:
%   phi        = coefficients of autoregressive model (with leading 1)
%   rho        = partial autocorrelations of generated model
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function [phi, rho] = ar_GenerateUniformCoefficients(p)

%% Generate partial autocorrelations
j = 1:p;
a = floor((j+1)/2);
b = floor(j/2)+1;
rho = 2*betarnd(a,b) - 1;

%% Transform to coefficients
phi = ar_PAC2Coef(rho);

return;