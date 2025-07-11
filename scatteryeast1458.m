clear;
% % 加载之前保存的数组
% load('yeast1458.mat');
% %'d', 'h','p','MDD_core','kcore','neighbor_Core','GAC','average_i'
% % 现在可以使用数组x和y了
% %% 创建散点图
% figure; % 创建一个新的图形窗口
% scatter(d,average_i, 'r+'); 
% xlabel('Degree');
% ylabel('Φ', 'Rotation', 0);
% 
% figure; % 创建一个新的图形窗口
% scatter(h,average_i, 'r+'); 
% xlabel('H-index');
% ylabel('Φ', 'Rotation', 0);
% 
% figure; % 创建一个新的图形窗口
% scatter(p,average_i, 'r+'); 
% xlabel('Pagerank');
% ylabel('Φ', 'Rotation', 0);
% 
% figure; % 创建一个新的图形窗口
% scatter(MDD_core,average_i, 'r+'); 
% xlabel('MDD');
% ylabel('Φ', 'Rotation', 0);
% 
% figure; % 创建一个新的图形窗口
% scatter(kcore,average_i, 'r+'); 
% xlabel('K-shell');
% ylabel('Φ', 'Rotation', 0);
% 
% figure; % 创建一个新的图形窗口
% scatter(neighbor_Core,average_i, 'r+'); 
% xlabel('NCC');
% ylabel('Φ', 'Rotation', 0);
% 
% figure; % 创建一个新的图形窗口
% scatter(GAC,average_i, 'r+'); 
% xlabel('GAC');
% ylabel('Φ', 'Rotation', 0);
% 
% figure; % 创建一个新的图形窗口
% scatter(KSGC,average_i, 'r+'); 
% xlabel('KSGC');
% ylabel('Φ', 'Rotation', 0);
load('yeast1458_sir1_30_Kendall.mat');
%'a00_degree', 'a00_h_index','a00_pagerank','a00_MDD_core','a00_GAC','a00_kcore','a00_neighbor_core','a00_KSGC'

x = 0.01:0.01:0.30;  % x轴数据

%figure; % 打开一个新的图形窗口
% 调整整个图的尺寸
fig = gcf;
fig.Position = [100, 100, 1000, 600]; % [left, bottom, width, height]
% 为每种算法绘制曲线，并设置颜色和标记
% 为每种算法绘制曲线，并设置颜色和标记

plot(x, a00_GAC, 's-', 'Color', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'DisplayName', 'GAC', 'MarkerSize', 6); % 黑色
hold on; 
plot(x, a00_h_index, 'o-', 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'H', 'MarkerSize', 6); % 红色
plot(x, a00_pagerank, '^-', 'Color', [0, 0, 1], 'MarkerFaceColor', [0, 0, 1], 'DisplayName', 'PRC', 'MarkerSize', 6); % 蓝色
plot(x, a00_MDD_core, 'v-', 'Color', [1, 0, 1], 'MarkerFaceColor', [1, 0, 1], 'DisplayName', 'MDD', 'MarkerSize', 6); % 暗洋红色
plot(x, a00_neighbor_core, 'd-', 'Color', [0, 128/255, 0], 'MarkerFaceColor',[0, 128/255, 0], 'DisplayName', 'NCC', 'MarkerSize', 6); % 暗黄色
plot(x, a00_degree, '<-', 'Color', [0, 0, 128/255], 'MarkerFaceColor', [0, 0, 128/255], 'DisplayName', 'DC', 'MarkerSize', 6); % 暗红色
plot(x, a00_kcore, '>-', 'Color', [128/255, 0, 1], 'MarkerFaceColor', [128/255, 0, 1], 'DisplayName', 'KS', 'MarkerSize', 6); % 深灰色
plot(x, a00_KSGC, 'h-', 'Color', [128/255, 0, 128/255], 'MarkerFaceColor', [128/255, 0, 128/255], 'DisplayName', 'KSGC', 'MarkerSize', 6); % 深粉红色

% 添加 a00_EKC_values 曲线（青绿色）
plot(x, a00_EKC_values, 'p-', 'Color', [0, 1, 1], 'MarkerFaceColor', [0, 1, 1], 'DisplayName', 'EKC', 'MarkerSize', 6); % 青色

% 添加 a00_neighbor_EKC_values 曲线（橙色）
plot(x, a00_neighbor_EKC_values, 'x-', 'Color', [1, 165/255, 0], 'MarkerFaceColor', [1, 165/255, 0], 'DisplayName', 'Neighbor EKC', 'MarkerSize', 6); % 橙色

% 添加 a00_extend_neighbor_EKC 曲线（紫蓝色）
plot(x, a00_extend_neighbor_EKC, '*-', 'Color', [75/255, 0, 130/255], 'MarkerFaceColor', [75/255, 0, 130/255], 'DisplayName', 'Extended Neighbor EKC', 'MarkerSize', 6); % 紫蓝色
% 添加图例
legend('show', 'Location', 'eastoutside');
%legend('show', 'Location', 'southwest');


% 添加标题和轴标签
title('yeast');
xlabel('\beta');
ylabel('\tau', 'Rotation', 0);

% 确保所有线都显示在图上
hold off;
