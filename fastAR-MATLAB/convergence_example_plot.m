% Load global figure settings
run('figure_settings.m');

%% Set up
close all
rng(2086)

ar_TestConvergence % Run the tests

%% Density 1
figure(1); clf; box on; hold on; grid on;
fig = gcf;
fig.Units = 'inches';
fig.Position = [5, 4, SMALL_FIG_WIDTH*1.25, SMALL_FIG_HEIGHT*1.1];
plot(log(n), mean(log(ME_sqrt_n)), mybluestyle, 'LineWidth', 1.1, 'color', [0,0,0]);
plot(log(n), mean(log(ME_2_sqrt_n)), '--'   , 'LineWidth', 1.1, 'color', [0.2,0.2,0.2]);
xlabel('$\log(n)$','FontSize',SMALL_FIG_AXIS_FONTSIZE,'Interpreter','Latex');
ylabel('$\log($ME$)$','FontSize',SMALL_FIG_AXIS_FONTSIZE,'Interpreter','Latex');
set(gca, 'FontSize', SMALL_FIG_AXIS_FONTSIZE);
set(gca, 'FontName', SMALL_FIG_FONTNAME);
set(gcf, 'color', 'w');
xlim([4.6,10.31]);

lgd = legend('~$p_n=\sqrt{n}$','~$p_n=2 \sqrt{n}$');
set(lgd, 'interpreter','latex','fontsize',LARGE_FIG_LEGEND_FONTSIZE);