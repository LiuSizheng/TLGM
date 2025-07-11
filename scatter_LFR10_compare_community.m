clear;


load('LFR2000-k10_average_sir1_15.mat');
load('LFR2000-10_GAC_Louvain_LP.mat')
%'a00_degree', 'a00_h_index','a00_pagerank','a00_MDD_core','a00_GAC','a00_kcore','a00_neighbor_core','a00_KSGC','a00_GRAD','a00_GRACC','a00_LRAD','a00_LRACC'
for i=1:15
    a00_GAC(i)=corr(average_iall{i}',GAC','type','Kendall');
    a00_GAC_Louvain(i)=corr(average_iall{i}',GAC_Louvain','type','Kendall');
    a00_GAC_LP(i)=corr(average_iall{i}',GAC_LP','type','Kendall');
end
x = 0.01:0.01:0.15;  % x轴数据

figure; % 打开一个新的图形窗口


plot(x, a00_GAC, 's-', 'Color', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'DisplayName', 'GAC', 'MarkerSize', 6); % 黑色
hold on; 
plot(x, a00_GAC_Louvain, 'o-', 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'GAC-Louvain', 'MarkerSize', 6); % 红色
plot(x, a00_GAC_LP, '^-', 'Color', [0, 0, 1], 'MarkerFaceColor', [0, 0, 1], 'DisplayName', 'GAC-LP', 'MarkerSize', 6); % 蓝色

% 添加图例
% legend('show', 'Location', 'eastoutside');
% legend('show', 'Location', 'southwest');


% 添加标题和轴标签
title('LFR,k=10');
xlabel('\beta');
ylabel('\tau', 'Rotation', 0);

% 确保所有线都显示在图上
hold off;
