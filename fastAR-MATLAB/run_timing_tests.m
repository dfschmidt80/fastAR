%% Run all the timing tests to recreate the table in
%
%   "Extremely Fast Maximum Likelihood Estimation of High-Order Autoregressive Models"
%   D. F. Schmidt and E. Makalic, Journal of Time Series Analysis (to appear)
%
% Note: the absolute timings will obviously differ depending on machine
% type, machine load, etc.

[s_005_a, s_0005_a, s_0001_a, delta_L_a] = ar_FastARTest(5, 10, 50, 1e2);
[s_005_b, s_0005_b, s_0001_b, delta_L_b] = ar_FastARTest(15, 20, 100, 1e2);
[s_005_c, s_0005_c, s_0001_c, delta_L_c] = ar_FastARTest(35, 50, 200, 1e2);
[s_005_d, s_0005_d, s_0001_d, delta_L_d] = ar_FastARTest(80, 100, 50, 1e2);