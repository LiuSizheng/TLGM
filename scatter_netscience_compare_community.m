clear;


load('netscience_GAC_GRA_sir1_26_Kendall.mat');
load('netscience_GAC_Louvain_LP.mat')
%'a00_degree', 'a00_h_index','a00_pagerank','a00_MDD_core','a00_GAC','a00_kcore','a00_neighbor_core','a00_KSGC','a00_GRAD','a00_GRACC','a00_LRAD','a00_LRACC'
for i=1:26
    a00_GAC_Louvain(i)=corr(average_iall{i}',GAC_Louvain','type','Kendall');
    a00_GAC_LP(i)=corr(average_iall{i}',GAC_LP','type','Kendall');
end
x = 0.01:0.01:0.26;  % x轴数据

figure; % 打开一个新的图形窗口

% 为每种算法绘制曲线，并设置颜色和标记
% 为每种算法绘制曲线，并设置颜色和标记

plot(x, a00_GAC, 's-', 'Color', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'DisplayName', 'GAC', 'MarkerSize', 6); % 黑色
hold on; 
plot(x, a00_GAC_Louvain, 'o-', 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'GAC-Louvain', 'MarkerSize', 6); % 红色
plot(x, a00_GAC_LP, '^-', 'Color', [0, 0, 1], 'MarkerFaceColor', [0, 0, 1], 'DisplayName', 'GAC-LP', 'MarkerSize', 6); % 蓝色
% plot(x, a00_h_index, 'o-', 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'H', 'MarkerSize', 6); % 红色
% plot(x, a00_pagerank, '^-', 'Color', [0, 0, 1], 'MarkerFaceColor', [0, 0, 1], 'DisplayName', 'PRC', 'MarkerSize', 6); % 蓝色
% plot(x, a00_MDD_core, 'v-', 'Color', [1, 0, 1], 'MarkerFaceColor', [1, 0, 1], 'DisplayName', 'MDD', 'MarkerSize', 6); % 暗洋红色
% plot(x, a00_neighbor_core, 'd-', 'Color', [0, 128/255, 0], 'MarkerFaceColor',[0, 128/255, 0], 'DisplayName', 'NCC', 'MarkerSize', 6); % 暗黄色
% plot(x, a00_degree, '<-', 'Color', [0, 0, 128/255], 'MarkerFaceColor', [0, 0, 128/255], 'DisplayName', 'DC', 'MarkerSize', 6); % 暗红色
% plot(x, a00_kcore, '>-', 'Color', [128/255, 0, 1], 'MarkerFaceColor', [128/255, 0, 1], 'DisplayName', 'KS', 'MarkerSize', 6); % 深灰色
% plot(x, a00_KSGC, 'h-', 'Color', [128/255, 0, 128/255], 'MarkerFaceColor', [128/255, 0, 128/255], 'DisplayName', 'KSGC', 'MarkerSize', 6); % 深粉红色
% plot(x, a00_GRAD, 'p-', 'Color', [0, 0.5, 0], 'MarkerFaceColor', [0, 0.5, 0], 'DisplayName', 'GRAD', 'MarkerSize', 6); % 深绿色，五角星
% plot(x, a00_GRACC, '+-', 'Color', [1, 0.5, 0], 'MarkerFaceColor', [1, 0.5, 0], 'DisplayName', 'GRACC', 'MarkerSize', 6); % 橙色，加号
% plot(x, a00_LRAD, '*-', 'Color', [0.5, 0.5, 0], 'MarkerFaceColor', [0.5, 0.5, 0], 'DisplayName', 'LRAD', 'MarkerSize', 6); % 黄色，星形
% plot(x, a00_LRACC, 'x-', 'Color', [0, 0.75, 0.75], 'MarkerFaceColor', [0, 0.75, 0.75], 'DisplayName', 'LRACC', 'MarkerSize', 6); % 青色，叉号

% 添加图例
%legend('show', 'Location', 'eastoutside');
legend('show', 'Location', 'southeast');


% 添加标题和轴标签
title('netscience');
xlabel('\beta');
ylabel('\tau', 'Rotation', 0);

% 确保所有线都显示在图上
hold off;
