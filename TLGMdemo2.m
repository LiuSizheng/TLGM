clear;% 清除工作区变量
dataset_name = 'protain';  % 加载数据集并处理成邻接矩阵examplenetwork
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
%[D,aver_D]=Aver_Path_Length(mixedsig);
neighbor_Core=zeros(N, 1);
for i=1:N
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    for j=1:len1
        neighbor_Core(i)=core(B(j))+neighbor_Core(i);
    end
end
% 定义计算半全局三角中心性的函数
% 初始化三角形数量向量
triangle_count = zeros(N, 1);
% 计算每个节点的三角形数量
for i = 1:N
    % 获取节点 i 的邻居节点
    neighbors = find(mixedsig(i, :));
    num_neighbors = length(neighbors);
    for j = 1:num_neighbors
        for k = j+1:num_neighbors
            % 检查邻居节点之间是否有连接
            if mixedsig(neighbors(j), neighbors(k)) == 1
                triangle_count(i) = triangle_count(i) + 1;
            end
        end
    end
end

%triangle_count_d=(triangle_count)/max(triangle_count)+d'/max(d');
triangle_count_d=(triangle_count)+d';
% 初始化半全局三角中心性向量
stc = zeros(N, 1);
for i = 1:N
    stc(i)=triangle_count_d(i);
    first_neighbors = find(mixedsig(i, :));
    triangle_count_d_1st = triangle_count_d(first_neighbors);
    stc(i)=stc(i)+sum(triangle_count_d_1st)/2;
    second_neighbors = [];
    for j = first_neighbors
        second_neighbors = [second_neighbors, find(mixedsig(j, :))];
    end
    second_neighbors = unique(second_neighbors);
    % 排除节点 i 本身和一阶邻居
    second_neighbors = setdiff(second_neighbors, [i, first_neighbors]);
    triangle_count_d_2nd = triangle_count_d(second_neighbors);
    stc(i)=stc(i)+sum(triangle_count_d_2nd)/4;
end

estc = zeros(N, 1);
for i = 1:N
    % 获取节点 i 的邻居节点
    neighbors = find(mixedsig(i, :));
    for j = 1:length(neighbors)
        % 计算半全局三角中心性
        estc(i) = estc(i) + stc(neighbors(j));
    end
end
% 计算新的中心性度量
%r = 3;% 设置距离阈值为3
k_stc=estc.*neighbor_Core;

% -------------------------------------------------------------------------
% 优化部分：使用邻接矩阵的幂来计算路径数量
% -------------------------------------------------------------------------

% 计算邻接矩阵的幂，用于路径计数
sparse_A = sparse(mixedsig);  % 确保使用稀疏矩阵以提高效率
A_squared = sparse_A^2;       % 二阶路径数量
A_cubed = sparse_A^3;         % 三阶路径数量

% 预计算节点的一阶、二阶、三阶邻居以提高效率
first_neighbors_cell = cell(N, 1);
second_neighbors_cell = cell(N, 1);
third_neighbors_cell = cell(N, 1);

for i = 1:N
    % 一阶邻居
    first_neighbors_cell{i} = find(mixedsig(i, :));
    
    % 二阶邻居 - 使用A^2并排除一阶邻居和节点自身
    potential_second = find(A_squared(i, :) > 0);
    second_neighbors_cell{i} = setdiff(potential_second, [i, first_neighbors_cell{i}]);
    
    % 三阶邻居 - 使用A^3并排除一阶、二阶邻居和节点自身
    potential_third = find(A_cubed(i, :) > 0);
    third_neighbors_cell{i} = setdiff(potential_third, [i, first_neighbors_cell{i}, second_neighbors_cell{i}]);
end

% 初始化DC_stc向量
DC_stc = zeros(N, 1);
% 初始化全局缓存矩阵，存储已计算的节点对电阻距离
R_global = NaN(N, N);  % 初值为NaN，表示未计算

for i = 1:N
    % 获取节点i的一阶、二阶、三阶邻居
    first_neighbors = first_neighbors_cell{i};
    second_neighbors = second_neighbors_cell{i};
    third_neighbors = third_neighbors_cell{i};
    
    % 合并所有邻居
    all_neighbors = union(union(first_neighbors, second_neighbors), third_neighbors);
    resistace_distance = zeros(N, 1);
    
    for r = all_neighbors
        if isnan(R_global(i, r))
            % 使用邻接矩阵的幂计算节点i和r之间的路径数量
            paths_1 = mixedsig(i, r);              % 长度为1的路径数量
            paths_2 = A_squared(i, r);             % 长度为2的路径数量
            paths_3 = A_cubed(i, r);               % 长度为3的路径数量
            
            % 计算等效电阻
            resist_sum = paths_1/1 + paths_2/2 + paths_3/3;
            
            % 更新缓存矩阵（利用对称性）
            temp_resist = 1 / resist_sum;
            R_global(i, r) = temp_resist;
            R_global(r, i) = temp_resist;
        end
        resistace_distance(r) = R_global(i, r);
    end
    
    % 计算与节点i距离小于等于3的节点的DC_stc
    DC_stc_sum = 0;
    % 一阶邻居
    for j1 = first_neighbors
        DC_stc_sum = DC_stc_sum + (k_stc(i) * k_stc(j1)) / resistace_distance(j1)^2;
    end
    % 二阶邻居
    for j2 = second_neighbors
        DC_stc_sum = DC_stc_sum + (k_stc(i) * k_stc(j2)) / resistace_distance(j2)^2;
    end
    % 三阶邻居
    for j3 = third_neighbors
        DC_stc_sum = DC_stc_sum + (k_stc(i) * k_stc(j3)) / resistace_distance(j3)^2;
    end
    % 将计算结果赋给DC_stc向量
    DC_stc(i) = DC_stc_sum;
end

% 打印前10个节点的TLGM值以验证算法
disp('前10个节点的TLGM值:');
disp(DC_stc(1:10));