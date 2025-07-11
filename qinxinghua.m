clear;% 清除工作区变量
dataset_name = 'email';  % 加载数据集并处理成邻接矩阵examplenetwork
A = load([dataset_name, '.txt']);
% facebook324 infectious protain CA-GrQc 1LFR2000-k5 
% LFR2000-10 LFR2000-k15 netscience 0.06-0.20 yeast1458 0.07-0.21
TT = A(:, 1:2); % 取数据集的前两列
mixedsig = zeros(max(max(TT)));
N = size(mixedsig, 1); % 获取邻接矩阵的节点数
len = length(TT); 
for i = 1:len
    mixedsig(TT(i, 1), TT(i, 2)) = 1;
    mixedsig(TT(i, 2), TT(i, 1)) = 1;
end
kk = sparse(mixedsig);
[a, b] = components(kk);
[B] = largestcomponent(mixedsig);
mixedsig = mixedsig(B, B);% 只保留最大连通分量对应的邻接矩阵
N = size(mixedsig, 1); % 获取邻接矩阵的节点数
A=mixedsig;
b=sum(A,2);
k1=sum(b)/N;
k2=sum(b'*b)/N;
irate=k1/k2;%传播率阈值
ccfs = sum(clustering_coefficients(sparse(mixedsig)))/N; %网络聚类系数
dmean=sum(degreea( mixedsig))/N; %网络平均度
d = degreea( mixedsig );
core=core_numbers(sparse(mixedsig));
[D,aver_D]=Aver_Path_Length(mixedsig);
neighbor_Core=zeros(N, 1);
for i=1:N
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    for j=1:len1
        neighbor_Core(i)=core(B(j))+neighbor_Core(i);
    end
end
L_pinv=calculate_laplacian_pinv(mixedsig);

R = zeros(N, N);     % 初始化 R 矩阵

for i = 1:N
    for j = 1:N
        % 按公式计算每个元素
        R(i, j) = abs(L_pinv(i, i) + L_pinv(j, j) - 2 * L_pinv(i, j));
    end
end
G_values = kshell_entropy_gravity_run_proposed(mixedsig, core, N,R);
matData = load('infectious_sir1_15_Kendall.mat');
average_iall = matData.average_iall;
for i=1:15
    a00_G_values(i)=corr(average_iall{i}',G_values,'type','Kendall');  
    % a00_DC(i)=corr(average_iall{i}',DC,'type','Kendall');  
end
x = 0.01:0.01:0.15;  % x轴数据
% 调整整个图的尺寸
fig = gcf;
fig.Position = [100, 100, 1000, 600]; % [left, bottom, width, height]
% 为每种算法绘制曲线，并设置颜色和标记
plot(x, matData.a00_GAC, 's-', 'Color', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'DisplayName', 'GLC', 'MarkerSize', 6); % 黑色
hold on;
plot(x, matData.a00_h_index, 'o-', 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'H', 'MarkerSize', 6); % 红色
plot(x, matData.a00_pagerank, '^-', 'Color', [0, 0, 1], 'MarkerFaceColor', [0, 0, 1], 'DisplayName', 'PRC', 'MarkerSize', 6); % 蓝色
plot(x, matData.a00_MDD_core, 'v-', 'Color', [1, 0, 1], 'MarkerFaceColor', [1, 0, 1], 'DisplayName', 'MDD', 'MarkerSize', 6); % 暗洋红色
plot(x, matData.a00_neighbor_core, 'd-', 'Color', [0, 128/255, 0], 'MarkerFaceColor',[0, 128/255, 0], 'DisplayName', 'NCC', 'MarkerSize', 6); % 暗黄色
plot(x, matData.a00_degree, '<-', 'Color', [0, 0, 128/255], 'MarkerFaceColor', [0, 0, 128/255], 'DisplayName', 'DC', 'MarkerSize', 6); % 暗红色
plot(x, matData.a00_kcore, '>-', 'Color', [128/255, 0, 1], 'MarkerFaceColor', [128/255, 0, 1], 'DisplayName', 'KS', 'MarkerSize', 6); % 深灰色
plot(x, matData.a00_KSGC, 'h-', 'Color', [128/255, 0, 128/255], 'MarkerFaceColor', [128/255, 0, 128/255], 'DisplayName', 'KSGC', 'MarkerSize', 6); % 深粉红色
% 添加 a00_DC_stc 曲线（青绿色）
plot(x, a00_G_values, 'p-', 'Color', [0, 1, 1], 'MarkerFaceColor', [0, 1, 1], 'DisplayName', 'DC_stc', 'MarkerSize', 6); % 青色
% 添加 a00_DC 曲线（紫蓝色）
% plot(x, a00_DC, 'x-', 'Color', [1, 165/255, 0], 'MarkerFaceColor', [1, 165/255, 0],  'DisplayName', 'density centrality', 'MarkerSize', 6); % 
ylim_vals = ylim;  % 获取当前 y 轴范围
plot([irate, irate], ylim_vals, 'k--', 'LineWidth', 1.5);
% 添加图例
%legend('show', 'Location', 'eastoutside');
legend('show', 'Location', 'eastoutside');
% 添加标题和轴标签
title(dataset_name);
xlabel('\beta');
ylabel('\tau', 'Rotation', 0);
% 确保所有线都显示在图上
hold off;