load('lambda_test.mat');%, 'average_kendall_all','GAC_iall'
%sir_database={'facebook324_average_sir1_20.mat','netscience_average_sir1_26.mat','infectious_average_sir1_20.mat','yeast1458_average_sir1_30.mat','protain_average_sir1_15.mat','CA-GrQc_average_sir1_15.mat'};
name={'facebook','netscience','infectious','yeast','protain','CA-GrQc'};
p=6;
figure;
x=0.05:0.05:1;
plot(x, cell2mat(average_kendall_all{p}), '^-', 'Color', [0, 0, 1], 'MarkerFaceColor', [0, 0, 1], 'DisplayName', 'PRC', 'MarkerSize', 6); % 黑色
% 添加标题和轴标签
title(name{p});
xlabel('$\lambda$', 'Interpreter', 'latex');
ylabel('average kendall''s $\textless\tau\textgreater$', 'Interpreter', 'latex');
%yticks([0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92]); % 设置y轴的刻度
ylim([0.77 0.79]); % 设置y轴的显示范围
