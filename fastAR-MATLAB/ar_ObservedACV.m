%% V1.0
%% Efficiently compute the observed autocovariances
% 
% function acf = ar_ObservedACV(y, nlags)
%
% Parameters:
%   y          = the time series
%   nlags      = the number of lags to compute
%
% Returns:
%   acf        = observed autocovariances
%
% Copyright (C) Daniel F Schmidt and Enes Makalic
%
function acf = ar_ObservedACV(y, nlags)

nFFT = 2^(nextpow2(length(y))+1);
F = fft(y,nFFT);
F = F.*conj(F);
acf = ifft(F);
acf = acf(1:(nlags+1)); % Retain nonnegative lags
acf = real(acf);

acf = acf/length(y);

end