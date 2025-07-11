% 清除工作区变量
clear;
% 加载数据集并处理成邻接矩阵
% 定义数据集名称
dataset_name = 'CA-GrQc';
A = load([dataset_name, '.txt']);%yeast1458 facebook324 protain netscience infectious CA-GrQc 1LFR2000-k5 LFR2000-10 LFR2000-k15
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

% 计算节点的 k-壳值
core = core_numbers(sparse(mixedsig)); 

first_neighbor_k=zeros(N,1);
second_neighbor_k=zeros(N,1);
for i = 1:N
    % 计算一阶邻居的 k-壳熵值 EK1
    first_neighbors = find(mixedsig(i, :)); 

    k_shell_values_1st = core(first_neighbors); 
    first_neighbor_k(i)=length(unique(k_shell_values_1st));
        
    second_neighbors = [];
    for j = first_neighbors
        second_neighbors = [second_neighbors, find(mixedsig(j, :))]; 
    end
    second_neighbors = unique(second_neighbors); 
    % 排除节点 i 本身和一阶邻居
    second_neighbors = setdiff(second_neighbors, [i, first_neighbors]);
    k_shell_values_2nd = core(second_neighbors); 
    second_neighbor_k(i)=length(unique(k_shell_values_2nd));

end

% 绘图部分
% 计算整个网络的最大 k-壳值作为横坐标范围
max_k = max(core); % 使用 core 的最大值

% 初始化统计数组，从 1 到 max_k（无孤立节点）
k_range = 1:max_k; % 从 1 开始，因为最大连通分量中没有孤立节点
first_count = zeros(size(k_range)); % 一阶邻居的计数
second_count = zeros(size(k_range)); % 二阶邻居的计数

% 统计每个 k 值的节点个数
for k = k_range
    first_count(k) = sum(first_neighbor_k == k); % k=1 对应索引 1
    second_count(k) = sum(second_neighbor_k == k);
end

% 绘制图形
figure; % 创建新图形窗口
plot(k_range, first_count, 'b-o', 'LineWidth', 2, 'DisplayName', 'First-order neighbors'); % 蓝色线表示一阶邻居
hold on;
plot(k_range, second_count, 'r-o', 'LineWidth', 2, 'DisplayName', 'Second-order neighbors'); % 红色线表示二阶邻居
hold off;

% 设置图形属性
xlabel('k-shell value (k)');
ylabel('Number of nodes');
title(dataset_name);
legend('show');
grid on;