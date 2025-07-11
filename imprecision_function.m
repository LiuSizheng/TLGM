clear;% 清除工作区变量
dataset_name = 'power-US-Grid';  % 加载数据集并处理成邻接矩阵examplenetwork
A = load([dataset_name, '.txt']);
% facebook324 USAir97 email infectious protain CA-GrQc 
% power-US-Grid   1LFR2000-k5 
% LFR2000-10 LFR2000-k15 netscience    yeast1458
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
% neighbor_Core=zeros(N, 1);
% for i=1:N
%     B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
%     len1=length(B);%len1表示数组B的长度
%     for j=1:len1
%         neighbor_Core(i)=core(B(j))+neighbor_Core(i);
%     end
% end
% % 定义计算半全局三角中心性的函数
% % 初始化三角形数量向量
% triangle_count = zeros(N, 1);
% % 计算每个节点的三角形数量
% for i = 1:N
%     % 获取节点 i 的邻居节点
%     neighbors = find(mixedsig(i, :));
%     num_neighbors = length(neighbors);
%     for j = 1:num_neighbors
%         for k = j+1:num_neighbors
%             % 检查邻居节点之间是否有连接
%             if mixedsig(neighbors(j), neighbors(k)) == 1
%                 triangle_count(i) = triangle_count(i) + 1;
%             end
%         end
%     end
% end
% 
% %triangle_count_d=(triangle_count)/max(triangle_count)+d'/max(d');
% triangle_count_d=(triangle_count)+d';
% % 初始化半全局三角中心性向量
% stc = zeros(N, 1);
% for i = 1:N
%     stc(i)=triangle_count_d(i);
%     first_neighbors = find(mixedsig(i, :));
%     triangle_count_d_1st = triangle_count_d(first_neighbors);
%     stc(i)=stc(i)+sum(triangle_count_d_1st)/2;
%     second_neighbors = [];
%     for j = first_neighbors
%         second_neighbors = [second_neighbors, find(mixedsig(j, :))];
%     end
%     second_neighbors = unique(second_neighbors);
%     % 排除节点 i 本身和一阶邻居
%     second_neighbors = setdiff(second_neighbors, [i, first_neighbors]);
%     triangle_count_d_2nd = triangle_count_d(second_neighbors);
%     stc(i)=stc(i)+sum(triangle_count_d_2nd)/4;
% end
% 
% estc = zeros(N, 1);
% for i = 1:N
%     % 获取节点 i 的邻居节点
%     neighbors = find(mixedsig(i, :));
%     for j = 1:length(neighbors)
%         % 计算半全局三角中心性
%         estc(i) = estc(i) + stc(neighbors(j));
%     end
% end
% % 计算新的中心性度量
% r = 3;% 设置距离阈值为3
% k_stc=estc.*neighbor_Core;
% % 初始化DC_stc向量
% DC_stc = zeros(N, 1);
% % 初始化全局缓存矩阵，存储已计算的节点对电阻距离，利用对称性降低计算量
% R_global = NaN(N, N);  % 初值为NaN，表示未计算
% %all_resitace_distance = cell(N, 1);
% %all_resitace_distance =zeros(N,N);
% for i = 1:N
%     % 获取节点i的一阶、二阶、三阶邻居
%     first_neighbors = find(mixedsig(i, :));  % 一阶邻居
%     second_neighbors = [];  % 二阶邻居
%     third_neighbors = [];   % 三阶邻居  
%     % 获取二阶邻居
%     for j1 = first_neighbors
%         second_neighbors = [second_neighbors, find(mixedsig(j1, :))]; % j1的邻居是i的二阶邻居
%     end
%     second_neighbors = unique(second_neighbors);
%     second_neighbors = setdiff(second_neighbors, [first_neighbors, i]); % 排除一阶邻居和节点i  
%     % 获取三阶邻居
%     for j2 = second_neighbors
%         third_neighbors = [third_neighbors, find(mixedsig(j2, :))]; % j2的邻居是i的三阶邻居
%     end
%     third_neighbors = unique(third_neighbors);
%     third_neighbors = setdiff(third_neighbors, [first_neighbors, second_neighbors, i]); % 排除已经考虑的邻居 
%     all_neighbors = union(union(first_neighbors, second_neighbors), third_neighbors);
%     resistace_distance=zeros(N,1);
%     for r = all_neighbors
%         if isnan(R_global(i, r))
%             resist_sum = 0;
%             % 一阶路径 (直接连接)
%             if mixedsig(i,r) == 1
%                 resist_sum = resist_sum + 1/1;  % 1步路径电阻倒数
%             end
% 
%             % 二阶路径
%             common_neighbors_2 = intersect(first_neighbors, find(mixedsig(r,:) == 1));
%             num_paths_2 = length(common_neighbors_2);
%             resist_sum = resist_sum + num_paths_2 * (1/2);  % 2步路径电阻倒数
% 
%             % 三阶路径
%             num_paths_3 = 0;
%             for l=first_neighbors
%                 if l~=r
%                     common_neighbors_lr=intersect(find(mixedsig(l,:) == 1), find(mixedsig(r,:) == 1));
%                     common_neighbors_lr = setdiff(common_neighbors_lr, i);
%                     num_paths_lr= length(common_neighbors_lr);
%                     num_paths_3 = num_paths_3 + num_paths_lr;
%                 end
%             end
%             resist_sum = resist_sum + num_paths_3 * (1/3);  % 3步路径电阻倒数
%             %resistace_distance(r)=1/resist_sum;
%             % 更新缓存矩阵（利用对称性）
%             temp_resist = 1 / resist_sum;
%             R_global(i, r) = temp_resist;
%             R_global(r, i) = temp_resist;
%         end
%         resistace_distance(r) = R_global(i, r);
%         %all_resitace_distance(i, r) = resistace_distance(r);
%     end
%     %all_resitace_distance{i}=resistace_distance;
% 
%     % 计算与节点i距离小于等于3的节点的DC_stc
%     DC_stc_sum = 0;
%     % 一阶邻居
%     for j1 = first_neighbors
%         DC_stc_sum = DC_stc_sum + (k_stc(i) * k_stc(j1)) / resistace_distance(j1)^2; % 距离为1
%     end
%     % 二阶邻居
%     for j2 = second_neighbors
%         DC_stc_sum = DC_stc_sum + (k_stc(i) * k_stc(j2)) / resistace_distance(j2)^2; % 距离为2，计算时平方
%     end
%     % 三阶邻居
%     for j3 = third_neighbors
%         DC_stc_sum = DC_stc_sum + (k_stc(i) * k_stc(j3)) / resistace_distance(j3)^2; % 距离为3，计算时平方
%     end
%     % 将计算结果赋给DC_stc向量
%     DC_stc(i) = DC_stc_sum;
% end

load('power-US-Grid.mat');  %yeast1458 facebook324 protain netscience infectious CA-GrQc 1LFR2000-k5 LFR2000-10 LFR2000-k15

% 将所有中心度数据统一转换为列向量形式 
%GAC = GAC(:);
KSGC = KSGC(:);
MDD_core = MDD_core(:);
average_i = average_i(:);
d = d(:);  % 度中心度
h = h(:);  % H指数中心度
kcore = kcore(:);  % K核中心度
neighbor_Core = neighbor_Core(:);  % 邻域中心度
p = p(:);  % PageRank中心度
DC_stc=DC_stc(:);
% EKC_values=EKC_values(:);
% neighbor_EKC=neighbor_EKC(:);
% extend_neighbor_EKC=extend_neighbor_EKC(:);
% 获取中心度的名称和数据
%centrality_names = {'GAC', 'KSGC', 'MDD', 'Degree', 'H_index', 'K_core', 'Neighbor_Core', 'PageRank','EKC','Neighbor EKC','Extend Neighbor EKC'};
centrality_data = {d,p,h,kcore,MDD_core,neighbor_Core,KSGC,DC_stc};
num_centralities = length(centrality_data);

% 参数设置
num_nodes = length(average_i);  % 网络节点总数
top_p_nodes = round(0.2 * num_nodes);  % 前 20% 的节点数量

% 计算基准传播效率 M_eff
[~, sorted_indices] = sort(average_i, 'descend');
M_eff = cumsum(average_i(sorted_indices(1:top_p_nodes))) ./ (1:top_p_nodes)';  % 基准传播效率随着前 k 个节点的均值

% 初始化结果存储矩阵
imprecision_values = zeros(top_p_nodes, num_centralities);

% 逐个计算每种中心度在不同节点数量下的不确定函数
for i = 1:num_centralities
    centrality = centrality_data{i};  % 当前中心度数据
    
    % 根据当前中心度对节点进行排序，计算前 k 个节点的不确定函数
    [~, sorted_indices] = sort(centrality, 'descend');
    sorted_infections = average_i(sorted_indices);  % 当前中心度排序后的感染效果
    
    % 计算不确定函数曲线 M_p
    M_p = cumsum(sorted_infections(1:top_p_nodes)) ./ (1:top_p_nodes)';  % 传播效率 M_p(k)
    imprecision_values(:, i) = 1 - (M_p ./ M_eff);  % 不确定函数 ε(k)
end

% 绘制不确定函数变化曲线
% figure;
% x=1:(top_p_nodes);
% x=x/num_nodes;
% plot(x, imprecision_values(:, 1), 's-', 'Color', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'DisplayName', 'GAC', 'MarkerSize', 6); % 黑色
% hold on; 
% plot(x, imprecision_values(:, 5), 'o-', 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'H', 'MarkerSize', 6); % 红色
% plot(x, imprecision_values(:, 8), '^-', 'Color', [0, 0, 1], 'MarkerFaceColor', [0, 0, 1], 'DisplayName', 'PRC', 'MarkerSize', 6); % 蓝色
% plot(x, imprecision_values(:, 3), 'v-', 'Color', [1, 0, 1], 'MarkerFaceColor', [1, 0, 1], 'DisplayName', 'MDD', 'MarkerSize', 6); % 暗洋红色
% plot(x, imprecision_values(:, 7), 'd-', 'Color', [0, 128/255, 0], 'MarkerFaceColor',[0, 128/255, 0], 'DisplayName', 'NCC', 'MarkerSize', 6); % 暗黄色
% plot(x, imprecision_values(:, 4), '<-', 'Color', [0, 0, 128/255], 'MarkerFaceColor', [0, 0, 128/255], 'DisplayName', 'DC', 'MarkerSize', 6); % 暗红色
% plot(x, imprecision_values(:, 6), '>-', 'Color', [128/255, 0, 1], 'MarkerFaceColor', [128/255, 0, 1], 'DisplayName', 'KS', 'MarkerSize', 6); % 深灰色
% plot(x, imprecision_values(:, 2), 'h-', 'Color', [128/255, 0, 128/255], 'MarkerFaceColor', [128/255, 0, 128/255], 'DisplayName', 'KSGC', 'MarkerSize', 6); % 深粉红色

%%
% 生成20个均匀分布的采样点
num_points = 20;
x = linspace(1, top_p_nodes, num_points) / num_nodes;  % 在前20%节点内取20个点
imprecision_sampled = zeros(num_points, num_centralities);

for i = 1:num_centralities
    % 对每个中心度的不确定函数进行插值，生成20个点
    imprecision_sampled(:, i) = interp1(1:top_p_nodes, imprecision_values(:, i), linspace(1, top_p_nodes, num_points));
end

% % 绘制不确定函数变化曲线
% figure;
% plot(x, imprecision_sampled(:, 1), 's-', 'Color', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'DisplayName', 'GLC', 'MarkerSize', 6); % 黑色
% hold on; 
% plot(x, imprecision_sampled(:, 5), 'o-', 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'H', 'MarkerSize', 6); % 红色
% plot(x, imprecision_sampled(:, 8), '^-', 'Color', [0, 0, 1], 'MarkerFaceColor', [0, 0, 1], 'DisplayName', 'PRC', 'MarkerSize', 6); % 蓝色
% plot(x, imprecision_sampled(:, 3), 'v-', 'Color', [1, 0, 1], 'MarkerFaceColor', [1, 0, 1], 'DisplayName', 'MDD', 'MarkerSize', 6); % 暗洋红色
% plot(x, imprecision_sampled(:, 7), 'd-', 'Color', [0, 128/255, 0], 'MarkerFaceColor',[0, 128/255, 0], 'DisplayName', 'NCC', 'MarkerSize', 6); % 暗绿色
% plot(x, imprecision_sampled(:, 4), '<-', 'Color', [0, 0, 128/255], 'MarkerFaceColor', [0, 0, 128/255], 'DisplayName', 'DC', 'MarkerSize', 6); % 深蓝色
% plot(x, imprecision_sampled(:, 6), '>-', 'Color', [128/255, 0, 1], 'MarkerFaceColor', [128/255, 0, 1], 'DisplayName', 'KS', 'MarkerSize', 6); % 深紫色
% plot(x, imprecision_sampled(:, 2), 'h-', 'Color', [128/255, 0, 128/255], 'MarkerFaceColor', [128/255, 0, 128/255], 'DisplayName', 'KSGC', 'MarkerSize', 6); % 深粉红色
% plot(x, imprecision_sampled(:, 9), 'p-', 'Color', [0, 1, 1], 'MarkerFaceColor', [0, 1, 1], 'DisplayName', 'EKC', 'MarkerSize', 6);% 青色
% plot(x, imprecision_sampled(:, 10), 'x-', 'Color', [1, 165/255, 0], 'MarkerFaceColor', [1, 165/255, 0], 'DisplayName', 'Neighbor EKC', 'MarkerSize', 6);% 橙色
% plot(x, imprecision_sampled(:, 11), '*-', 'Color', [75/255, 0, 130/255], 'MarkerFaceColor', [75/255, 0, 130/255], 'DisplayName', 'Extended Neighbor EKC', 'MarkerSize', 6);% 紫蓝色



%%

% 添加图例
%legend('show', 'Location', 'eastoutside');
% legend('show', 'Location', 'northeast');
% 
% % 添加标题和轴标签
% title('LFR2000-k15');
% xlabel('$p$', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('$\varepsilon(p)$', 'Interpreter', 'latex', 'FontSize', 14);
% xlim([0 0.205]); % 限制 x 轴的范围
% % 确保所有线都显示在图上
% hold off;
% save('imprecision_function.mat','imprecision_sampled');
