clear;% 清除工作区变量
dataset_name = 'LFR2000-10';  % 加载数据集并处理成邻接矩阵examplenetwork
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
r = 3;% 设置距离阈值为3
k_stc=estc.*neighbor_Core;
% 初始化DC_stc向量
DC_stc = zeros(N, 1);
% 计算DC_stc
for i = 1:N
    % 获取节点i的一阶、二阶、三阶邻居
    first_neighbors = find(mixedsig(i, :));  % 一阶邻居
    second_neighbors = [];  % 二阶邻居
    third_neighbors = [];   % 三阶邻居  
    % 获取二阶邻居
    for j1 = first_neighbors
        second_neighbors = [second_neighbors, find(mixedsig(j1, :))]; % j1的邻居是i的二阶邻居
    end
    second_neighbors = unique(second_neighbors);
    second_neighbors = setdiff(second_neighbors, [first_neighbors, i]); % 排除一阶邻居和节点i  
    % 获取三阶邻居
    for j2 = second_neighbors
        third_neighbors = [third_neighbors, find(mixedsig(j2, :))]; % j2的邻居是i的三阶邻居
    end
    third_neighbors = unique(third_neighbors);
    third_neighbors = setdiff(third_neighbors, [first_neighbors, second_neighbors, i]); % 排除已经考虑的邻居 
    % 计算与节点i距离小于等于3的节点的DC_stc
    DC_stc_sum = 0;
    % 一阶邻居
    for j1 = first_neighbors
        DC_stc_sum = DC_stc_sum + (k_stc(i) * k_stc(j1)) / 1; % 距离为1
    end
    % 二阶邻居
    for j2 = second_neighbors
        DC_stc_sum = DC_stc_sum + (k_stc(i) * k_stc(j2)) / 4; % 距离为2，计算时平方
    end
    % 三阶邻居
    for j3 = third_neighbors
        DC_stc_sum = DC_stc_sum + (k_stc(i) * k_stc(j3)) / 9; % 距离为3，计算时平方
    end
    % 将计算结果赋给DC_stc向量
    DC_stc(i) = DC_stc_sum;
end
% M=zeros(1, 8);
% RankData=load('LFR2000-k15.mat');
% M(1)=calculate_monotonicity(RankData.d);
% M(2) = calculate_monotonicity(RankData.p);
% M(3)=calculate_monotonicity(RankData.h);
% M(4) = calculate_monotonicity(RankData.kcore);
% M(5)=calculate_monotonicity(RankData.MDD_core);
% M(6) = calculate_monotonicity(RankData.neighbor_Core);
% M(7)=calculate_monotonicity(RankData.KSGC);
% M(8) = calculate_monotonicity(DC_stc);
% ccdf=cell(8, 1);
% ccdf{1}=CCDF(RankData.d);
% ccdf{2}=CCDF(RankData.p);
% ccdf{3}=CCDF(RankData.h);
% ccdf{4}=CCDF(RankData.kcore);
% ccdf{5}=CCDF(RankData.MDD_core);
% ccdf{6}=CCDF(RankData.neighbor_Core);
% ccdf{7}=CCDF(RankData.KSGC);
% ccdf{8}=CCDF(DC_stc);

matData = load('LFR2000-10_sir1_15_Kendall.mat');
average_iall = matData.average_iall;
for i=1:15
    a00_DC_stc(i)=corr(average_iall{i}',DC_stc,'type','Kendall');  
    % a00_DC(i)=corr(average_iall{i}',DC,'type','Kendall');  
end
% x = 0.01:0.01:0.15;  % x轴数据
% % 调整整个图的尺寸
% fig = gcf;
% fig.Position = [100, 100, 1000, 600]; % [left, bottom, width, height]
% % 为每种算法绘制曲线，并设置颜色和标记
% plot(x, matData.a00_GAC, 's-', 'Color', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'DisplayName', 'GLC', 'MarkerSize', 6); % 黑色
% hold on;
% plot(x, matData.a00_h_index, 'o-', 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'H', 'MarkerSize', 6); % 红色
% plot(x, matData.a00_pagerank, '^-', 'Color', [0, 0, 1], 'MarkerFaceColor', [0, 0, 1], 'DisplayName', 'PRC', 'MarkerSize', 6); % 蓝色
% plot(x, matData.a00_MDD_core, 'v-', 'Color', [1, 0, 1], 'MarkerFaceColor', [1, 0, 1], 'DisplayName', 'MDD', 'MarkerSize', 6); % 暗洋红色
% plot(x, matData.a00_neighbor_core, 'd-', 'Color', [0, 128/255, 0], 'MarkerFaceColor',[0, 128/255, 0], 'DisplayName', 'NCC', 'MarkerSize', 6); % 暗黄色
% plot(x, matData.a00_degree, '<-', 'Color', [0, 0, 128/255], 'MarkerFaceColor', [0, 0, 128/255], 'DisplayName', 'DC', 'MarkerSize', 6); % 暗红色
% plot(x, matData.a00_kcore, '>-', 'Color', [128/255, 0, 1], 'MarkerFaceColor', [128/255, 0, 1], 'DisplayName', 'KS', 'MarkerSize', 6); % 深灰色
% plot(x, matData.a00_KSGC, 'h-', 'Color', [128/255, 0, 128/255], 'MarkerFaceColor', [128/255, 0, 128/255], 'DisplayName', 'KSGC', 'MarkerSize', 6); % 深粉红色
% % 添加 a00_DC_stc 曲线（青绿色）
% plot(x, a00_DC_stc, 'p-', 'Color', [0, 1, 1], 'MarkerFaceColor', [0, 1, 1], 'DisplayName', 'DC_stc', 'MarkerSize', 6); % 青色
% % 添加 a00_DC 曲线（紫蓝色）
% % plot(x, a00_DC, 'x-', 'Color', [1, 165/255, 0], 'MarkerFaceColor', [1, 165/255, 0],  'DisplayName', 'density centrality', 'MarkerSize', 6); % 
% ylim_vals = ylim;  % 获取当前 y 轴范围
% plot([irate, irate], ylim_vals, 'k--', 'LineWidth', 1.5);
% % 添加图例
% %legend('show', 'Location', 'eastoutside');
% legend('show', 'Location', 'eastoutside');
% % 添加标题和轴标签
% title(dataset_name);
% xlabel('\beta');
% ylabel('\tau', 'Rotation', 0);
% % 确保所有线都显示在图上
% hold off;

