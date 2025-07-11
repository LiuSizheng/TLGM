
% 清除工作区变量
clear;
% 加载数据集并处理成邻接矩阵
dataset_name = 'LFR2000-k15';
A = load([dataset_name, '.txt']);%yeast1458 facebook324 protain netscience infectious CA-GrQc 1LFR2000-k5 LFR2000-10 LFR2000-k15
% A = A + 1; % 可根据需要对 A 进行操作
r = 0.2; % 可能是一个参数，具体作用需根据上下文确定
TT = A(:, 1:2); % 取数据集的前两列，假设表示边的两个节点
% 初始化邻接矩阵，矩阵大小根据节点的最大编号确定
mixedsig = zeros(max(max(TT)));
N = size(mixedsig, 1); % 获取邻接矩阵的节点数
len = length(TT); % 边的数量
% 将边信息填入邻接矩阵，由于是无向图，对称填充
for i = 1:len
    mixedsig(TT(i, 1), TT(i, 2)) = 1;
    mixedsig(TT(i, 2), TT(i, 1)) = 1;
end
kk = sparse(mixedsig);
% 计算图的连通分量，components 函数未给出，需要根据具体实现替换
[a, b] = components(kk);
% 获取最大连通分量，largestcomponent 函数未给出，需要根据具体实现替换
[B] = largestcomponent(mixedsig);
% 只保留最大连通分量对应的邻接矩阵
mixedsig = mixedsig(B, B);
N = size(mixedsig, 1); % 获取邻接矩阵的节点数
A=mixedsig;
N=size(A,1);
b=sum(A,2);
k1=sum(b)/N;
k2=sum(b'*b)/N;
irate=k1/k2;
% 计算节点的 k-壳值
core = core_numbers(sparse(mixedsig));
% 存储每个节点的 EKC 值
EKC_values = zeros(N, 1);
for i = 1:N
    % 计算一阶邻居的 k-壳熵值 EK1
    first_neighbors = find(mixedsig(i, :));
    if ~isempty(first_neighbors)
        k_shell_values_1st = core(first_neighbors);
        K_sum1 = sum(k_shell_values_1st);
        EK1 = -sum((k_shell_values_1st / K_sum1).* log(k_shell_values_1st / K_sum1));
        % 计算二阶邻居的 k-壳熵值 EK2
        second_neighbors = [];
        for j = first_neighbors
            second_neighbors = [second_neighbors, find(mixedsig(j, :))];
        end
        second_neighbors = unique(second_neighbors);
        % 排除节点 i 本身和一阶邻居
        second_neighbors = setdiff(second_neighbors, [i, first_neighbors]);
        k_shell_values_2nd = core(second_neighbors);
        K_sum2 = sum(k_shell_values_2nd);
        EK2 = -sum((k_shell_values_2nd / K_sum2).* log(k_shell_values_2nd / K_sum2));
        % 计算 lambda
        lambda = 1 / (1 + exp(-(EK2 - EK1)));
        % 计算节点的 k-壳多样性 EKC
        EKC_values(i) = EK1 + lambda * EK2;
    else
        EKC_values(i) = 0; % 对于孤立节点，EKC 值设为 0
    end
end
d = degreea( mixedsig );
[r,l]=size(mixedsig);
neighbor_EKC=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_EKC(i)=0;
    for j=1:len1
        neighbor_EKC(i)=EKC_values(B(j))+neighbor_EKC(i);
    end
end
[r,l]=size(mixedsig);
extend_neighbor_EKC=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    extend_neighbor_EKC(i)=0;
    for j=1:len1
        extend_neighbor_EKC(i)=neighbor_EKC(B(j))+extend_neighbor_EKC(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EDC_values = zeros(N, 1);
ks_max=max(core);
for i = 1:N
    % 计算一阶邻居的 k-壳熵值 EK1
    first_neighbors = find(mixedsig(i, :));
    if ~isempty(first_neighbors)
        k_shell_values_1st = core(first_neighbors);
        K_sum1 = sum(k_shell_values_1st);
        ED11 = -sum((k_shell_values_1st/ks_max).*(k_shell_values_1st / K_sum1).* log(k_shell_values_1st / K_sum1));
        kn1= accumarray(k_shell_values_1st(:), 1, [ks_max, 1]);
        ED12=0;
        for m = 1:ks_max
            if kn1(m) > 0
                ED12 = ED12 - (1 / (ks_max - m + 1)) * (kn1(m) / b(i)) * log(kn1(m) / b(i));
            end
        end

        % 计算二阶邻居的 k-壳熵值 EK2
        second_neighbors = [];
        for j = first_neighbors
            second_neighbors = [second_neighbors, find(mixedsig(j, :))];
        end
        second_neighbors = unique(second_neighbors);
        % 排除节点 i 本身和一阶邻居
        second_neighbors = setdiff(second_neighbors, [i, first_neighbors]);
        k_shell_values_2nd = core(second_neighbors);
        K_sum2 = sum(k_shell_values_2nd);
        ED21 = -sum((k_shell_values_2nd/ks_max).*(k_shell_values_2nd / K_sum2).* log(k_shell_values_2nd / K_sum2));
        kn2= accumarray(k_shell_values_2nd(:), 1, [ks_max, 1]);
        ED22=0;
        for n = 1:ks_max
            if kn2(n) > 0
                ED22 = ED22 - (1 / (ks_max - n + 1)) * (kn2(n) / length(second_neighbors)) * log(kn2(n) / length(second_neighbors));
            end
        end
        ED1=ED11+ED12;
        ED2=ED21+ED22;
        % 计算 lambda
        lambda = 1 / (1 + exp(-(ED2 - ED1)));
        % 计算节点的 k-壳多样性 EKC
        EDC_values(i) = ED1 +  lambda*ED2;
    else
        EDC_values(i) = 0; % 对于孤立节点，EKC 值设为 0
    end
end
d = degreea( mixedsig );
[r,l]=size(mixedsig);
neighbor_EDC=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_EDC(i)=0;
    for j=1:len1
        neighbor_EDC(i)=EDC_values(B(j))+neighbor_EDC(i);
    end
end
[r,l]=size(mixedsig);
extend_neighbor_EDC=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    extend_neighbor_EDC(i)=0;
    for j=1:len1
        extend_neighbor_EDC(i)=neighbor_EDC(B(j))+extend_neighbor_EDC(i);
    end
end

matData = load('LFR2000-k15_sir1_15_Kendall.mat');
average_iall = matData.average_iall;
for i=1:15

    a00_extend_neighbor_EKC(i)=corr(average_iall{i}',extend_neighbor_EKC','type','Kendall');

    a00_extend_neighbor_EDC(i)=corr(average_iall{i}',extend_neighbor_EDC','type','Kendall');
    a00_EDC_values(i)=corr(average_iall{i}',EDC_values,'type','Kendall');
    a00_neighbor_EDC_values(i)=corr(average_iall{i}',neighbor_EDC','type','Kendall');
    
end

x = 0.01:0.01:0.15;  % x轴数据

%figure; % 打开一个新的图形窗口

% 调整整个图的尺寸
fig = gcf;
fig.Position = [100, 100, 1000, 600]; % [left, bottom, width, height]
% 为每种算法绘制曲线，并设置颜色和标记
% 为每种算法绘制曲线，并设置颜色和标记

plot(x, matData.a00_GAC, 's-', 'Color', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'DisplayName', 'GAC', 'MarkerSize', 6); % 黑色
hold on;
plot(x, matData.a00_h_index, 'o-', 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'H', 'MarkerSize', 6); % 红色
plot(x, matData.a00_pagerank, '^-', 'Color', [0, 0, 1], 'MarkerFaceColor', [0, 0, 1], 'DisplayName', 'PRC', 'MarkerSize', 6); % 蓝色
plot(x, matData.a00_MDD_core, 'v-', 'Color', [1, 0, 1], 'MarkerFaceColor', [1, 0, 1], 'DisplayName', 'MDD', 'MarkerSize', 6); % 暗洋红色
plot(x, matData.a00_neighbor_core, 'd-', 'Color', [0, 128/255, 0], 'MarkerFaceColor',[0, 128/255, 0], 'DisplayName', 'NCC', 'MarkerSize', 6); % 暗黄色
plot(x, matData.a00_degree, '<-', 'Color', [0, 0, 128/255], 'MarkerFaceColor', [0, 0, 128/255], 'DisplayName', 'DC', 'MarkerSize', 6); % 暗红色
plot(x, matData.a00_kcore, '>-', 'Color', [128/255, 0, 1], 'MarkerFaceColor', [128/255, 0, 1], 'DisplayName', 'KS', 'MarkerSize', 6); % 深灰色
plot(x, matData.a00_KSGC, 'h-', 'Color', [128/255, 0, 128/255], 'MarkerFaceColor', [128/255, 0, 128/255], 'DisplayName', 'KSGC', 'MarkerSize', 6); % 深粉红色

% % 添加 a00_extend_neighbor_EDC 曲线（青绿色）
% plot(x, a00_extend_neighbor_EDC, 'p-', 'Color', [0, 1, 1], 'MarkerFaceColor', [0, 1, 1], 'DisplayName', 'Extended Neighbor EDC', 'MarkerSize', 6); % 青色
% 
% 
% % 添加 a00_extend_neighbor_EKC 曲线（紫蓝色）
% plot(x, a00_extend_neighbor_EKC, '*-', 'Color', [75/255, 0, 130/255], 'MarkerFaceColor', [75/255, 0, 130/255], 'DisplayName', 'Extended Neighbor EKC', 'MarkerSize', 6); % 紫蓝色

% 添加 a00_EKC_values 曲线（青绿色）
plot(x, a00_EDC_values, 'p-', 'Color', [0, 1, 1], 'MarkerFaceColor', [0, 1, 1], 'DisplayName', 'EDC', 'MarkerSize', 6); % 青色

% 添加 a00_neighbor_EKC_values 曲线（橙色）
plot(x, a00_neighbor_EDC_values, 'x-', 'Color', [1, 165/255, 0], 'MarkerFaceColor', [1, 165/255, 0], 'DisplayName', 'Neighbor EDC', 'MarkerSize', 6); % 橙色

% 添加 a00_extend_neighbor_EKC 曲线（紫蓝色）
plot(x, a00_extend_neighbor_EDC, '*-', 'Color', [75/255, 0, 130/255], 'MarkerFaceColor', [75/255, 0, 130/255], 'DisplayName', 'Extended Neighbor EDC', 'MarkerSize', 6); % 紫蓝色


% 添加 a00_extend_neighbor_EKC 曲线（紫蓝色）
plot(x, a00_extend_neighbor_EKC, '+-', 'Color', [0.5, 0.5, 0.5], 'MarkerFaceColor',  [0.5, 0.5, 0.5], 'DisplayName', 'Extended Neighbor EKC', 'MarkerSize', 6); % 

% 添加图例
%legend('show', 'Location', 'eastoutside');
legend('show', 'Location', 'eastoutside');


% 添加标题和轴标签
title(dataset_name);
xlabel('\beta');
ylabel('\tau', 'Rotation', 0);

% 确保所有线都显示在图上
hold off;
