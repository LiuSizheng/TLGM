
% 清除工作区变量
clear;
% 加载数据集并处理成邻接矩阵
A = load('yeast1458.txt');
% facebook324 infectious protain CA-GrQc 1LFR2000-k5 
% LFR2000-10 LFR2000-k15 netscience 0.06-0.20 yeast1458 0.07-0.21
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
h = h_index(mixedsig);
N = size(mixedsig, 1); % 获取邻接矩阵的节点数
A=mixedsig;
N=size(A,1);
b=sum(A,2);
k1=sum(b)/N;
k2=sum(b'*b)/N;
irate=k1/k2;
% 计算节点的 k-壳值
core = core_numbers(sparse(mixedsig));
r=sum(core)/sum(h);
m=[];
for i=1:N
    m(i)=core(i)+r*h(i)
end
%[G, G_plus] = gravity_centrality(mixedsig, core');
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
        % 计算二阶邻居
        second_neighbors = [];
        for neighbor1 = first_neighbors
            temp_neighbors = find(mixedsig(neighbor1, :));
            % 排除已经是一阶邻居和节点自身
            valid_neighbors = setdiff(temp_neighbors, [first_neighbors, i]);
            second_neighbors = union(second_neighbors, valid_neighbors);
        end
        
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
KEKC = calculate_pc(mixedsig, m, EKC_values);

[r,l]=size(mixedsig);
neighbor_KEKC=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_KEKC(i)=0;
    for j=1:len1
        neighbor_KEKC(i)=KEKC(B(j))+neighbor_KEKC(i);
    end
end
[r,l]=size(mixedsig);
extend_neighbor_KEKC=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    extend_neighbor_KEKC(i)=0;
    for j=1:len1
        extend_neighbor_KEKC(i)=neighbor_KEKC(B(j))+extend_neighbor_KEKC(i);
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

[r,l]=size(mixedsig);
neighbor_Core=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_Core(i)=0;
    for j=1:len1
        neighbor_Core(i)=core(B(j))+neighbor_Core(i);
    end
end
[r,l]=size(mixedsig);
extend_neighbor_core=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    extend_neighbor_core(i)=0;
    for j=1:len1
        extend_neighbor_core(i)=neighbor_Core(B(j))+extend_neighbor_core(i);
    end
end

matData = load('yeast1458_sir1_30_Kendall.mat');
average_iall = matData.average_iall;
for i=7:21
    a00_PC(i)=corr(average_iall{i}',KEKC,'type','Kendall');
    a00_neighbor_pc(i)=corr(average_iall{i}',neighbor_KEKC','type','Kendall');
    a00_extend_neighbor_pc(i)=corr(average_iall{i}',extend_neighbor_KEKC','type','Kendall');
end
% average_i=sir_main(mixedsig,1000,0.05) ;
%a00_IC=corr(average_i',IC,'type','Kendall');
%a00_HNPR=corr(average_i',HNPR,'type','Kendall');
%a00_HNPRF=corr(average_i',HNPRF,'type','Kendall');
%a00_neighbor_SCC=corr(average_i',neighbor_IC','type','Kendall');
%a00_extend_neighbor_SCC=corr(average_i',extend_neighbor_IC','type','Kendall');
% a00_KEKC=corr(average_i',KEKC,'type','Kendall');
% a00_neighbor_KEKC=corr(average_i',neighbor_KEKC','type','Kendall');
% a00_extend_neighbor_KEKC=corr(average_i',extend_neighbor_KEKC','type','Kendall');
% a00_degree=corr(average_i',d','type','Kendall');
%a00_SCC=corr(average_i',SCC,'type','Kendall');
%a00_neighbor_IC=corr(average_i',neighbor_IC','type','Kendall');
%a00_extend_neighbor_IC=corr(average_i',extend_neighbor_IC','type','Kendall');

%EKC_values=EKC_values';
%a00_EKC_values=corr(average_i',EKC_values','type','Kendall');
%a00_core=corr(average_i',core,'type','Kendall');
%a00_neighbor_EKC_values=corr(average_i',neighbor_EKC','type','Kendall');
%a00_extend_neighbor_EKC=corr(average_i',extend_neighbor_EKC','type','Kendall');
% a00_neighbor_Core=corr(average_i',neighbor_Core','type','Kendall');
% a00_extend_neighbor_core=corr(average_i',extend_neighbor_core','type','Kendall');


