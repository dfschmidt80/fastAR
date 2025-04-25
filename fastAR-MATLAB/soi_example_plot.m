% Load global figure settings
run('figure_settings.m');

%% Set up
close all

% Load in the SOI data and divide into training and testing data
df = importdata('DailySOI1887_1989Base.csv');

n_test = 14;
n = length(df.data);
n_train = n - n_test;

y = df.data(1:n_train,end);
yf = df.data(n_train+1:end,end);

%% Fit the models to the first n_train observations using ML
pmax = ceil(2.5*sqrt(n_train));
tic;rv = ar_FitNested(y,pmax,'reltol',0.0001);toc

% Compute model selection scores
[phi_hat, sigma2_hat, mu_hat, L, rho_hat, score] = ar_Select(rv, 'AICc');
aicc = score{:,3};
aicc = aicc - min(aicc);

% Compute the forecasts and calculate the forecast error
y_hat_ml = ar_ForecastAR(y, phi_hat, mu_hat, n_test);
fprintf('Forecast error (ML): %f\n', mean((y_hat_ml - yf).^2));

%% Generate the plot 
figure(1); clf; box on; hold on; grid on;
fig = gcf;
fig.Units = 'inches';
fig.Position = [5, 4, SMALL_FIG_WIDTH*1.25, SMALL_FIG_HEIGHT*1.1];

plot(0:pmax, aicc, mybluestyle, 'LineWidth', 1.1, 'color', [0,0,0]);

xlabel('$p$','FontSize',SMALL_FIG_AXIS_FONTSIZE,'Interpreter','Latex');
ylabel('AIC$_c$','FontSize',SMALL_FIG_AXIS_FONTSIZE,'Interpreter','Latex');
set(gca, 'FontSize', SMALL_FIG_AXIS_FONTSIZE);
set(gca, 'FontName', SMALL_FIG_FONTNAME);
set(gcf, 'color', 'w');
xlim([0,212]);
ylim([0,150]);

%export_fig 'soi-example.pdf' % Export it (for the paper)