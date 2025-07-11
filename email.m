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
% 初始化DC_stc向量
DC_stc = zeros(N, 1);
% 初始化全局缓存矩阵，存储已计算的节点对电阻距离，利用对称性降低计算量
R_global = NaN(N, N);  % 初值为NaN，表示未计算
%all_resitace_distance = cell(N, 1);
%all_resitace_distance =zeros(N,N);
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
    all_neighbors = union(union(first_neighbors, second_neighbors), third_neighbors);
    resistace_distance=zeros(N,1);
    for r = all_neighbors
        if isnan(R_global(i, r))
            resist_sum = 0;
            % 一阶路径 (直接连接)
            if mixedsig(i,r) == 1
                resist_sum = resist_sum + 1/1;  % 1步路径电阻倒数
            end
            
            % 二阶路径
            common_neighbors_2 = intersect(first_neighbors, find(mixedsig(r,:) == 1));
            num_paths_2 = length(common_neighbors_2);
            resist_sum = resist_sum + num_paths_2 * (1/2);  % 2步路径电阻倒数
            
            % 三阶路径
            num_paths_3 = 0;
            for l=first_neighbors
                if l~=r
                    common_neighbors_lr=intersect(find(mixedsig(l,:) == 1), find(mixedsig(r,:) == 1));
                    common_neighbors_lr = setdiff(common_neighbors_lr, i);
                    num_paths_lr= length(common_neighbors_lr);
                    num_paths_3 = num_paths_3 + num_paths_lr;
                end
            end
            resist_sum = resist_sum + num_paths_3 * (1/3);  % 3步路径电阻倒数
            %resistace_distance(r)=1/resist_sum;
            % 更新缓存矩阵（利用对称性）
            temp_resist = 1 / resist_sum;
            R_global(i, r) = temp_resist;
            R_global(r, i) = temp_resist;
        end
        resistace_distance(r) = R_global(i, r);
        %all_resitace_distance(i, r) = resistace_distance(r);
    end
    %all_resitace_distance{i}=resistace_distance;
    
    % 计算与节点i距离小于等于3的节点的DC_stc
    DC_stc_sum = 0;
    % 一阶邻居
    for j1 = first_neighbors
        DC_stc_sum = DC_stc_sum + (k_stc(i) * k_stc(j1)) / resistace_distance(j1)^2; % 距离为1
    end
    % 二阶邻居
    for j2 = second_neighbors
        DC_stc_sum = DC_stc_sum + (k_stc(i) * k_stc(j2)) / resistace_distance(j2)^2; % 距离为2，计算时平方
    end
    % 三阶邻居
    for j3 = third_neighbors
        DC_stc_sum = DC_stc_sum + (k_stc(i) * k_stc(j3)) / resistace_distance(j3)^2; % 距离为3，计算时平方
    end
    % 将计算结果赋给DC_stc向量
    DC_stc(i) = DC_stc_sum;
end
adjMat=mixedsig;
d = degreea( adjMat ); %度中心性

h = h_index(adjMat); %H指数

p=pagerank(adjMat,N);%PageRank

MDD_core=MDD(adjMat);%混合度分解法（MDD）

kcore=core_numbers(sparse(adjMat))'; %k-壳分解

[r,l]=size(adjMat);  %邻域核度NCC
neighbor_Core=[];
for i=1:r
    B=find(adjMat(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_Core(i)=0;
    for j=1:len1
        neighbor_Core(i)=kcore(B(j))+neighbor_Core(i);
    end
end

cij=max(kcore)-min(kcore);  %KSGC算法
for i=1:N
    disp(i);
    j=1;%第一层
    nei0=[i];%第0层
    nei1=find(adjMat(i,:));%算出节点i的直接邻居
    par=[];
    nei1=find(adjMat(i,:));
    len1=length(nei1);
    par1=[];
    for j=1:len1
        par1(j)=exp((kcore(i)-kcore(nei1(j)))/cij)*d(i)*d(nei1(j));
    end
    
    par11=sum(par1);
    j=2;%第二层
    %     nei2=neighbor(mixedsig,nei1');
    % a=neighbor(mixedsig,nei1);
    a=[];
    for v=1:length(nei1)
        %     rr=mixedsig(nei1(v));
        k=find(adjMat(nei1(v),:));
        a=union(k,a);
    end
    c=setdiff(a,nei0);%得到节点i的二跳邻居
    nei2=setdiff(c,nei1);%得到节点i的二跳邻居
    
    par2=[];
    for k=1:length(nei2)
        par2(k)=exp((kcore(i)-kcore(nei2(k)))/cij)*d(i)*d(nei2(k))/2^2;
        %         par2(k)=BBB(i,nei2(k))*(1/k1)^2+BBA(i,nei2(k))*(1/k1)^3;
    end
    par22=sum(par2);
    
    j=3;%第3层
    a=[];
    for v=1:length(nei2)
        %     rr=mixedsig(nei1(v));
        k=find(adjMat(nei2(v),:));
        a=union(k,a);
    end
    c=setdiff(a,union(nei0,nei1));%得到节点i的3跳邻居
    nei3=setdiff(c,nei2);
    % BBA4=mixedsig^4;
    par3=[];
    for k=1:length(nei3)
        par3(k)=exp((kcore(i)-kcore(nei3(k)))/cij)*d(i)*d(nei3(k))/3^2;
    end
    par33=sum(par3);
%          runim(i)=par11+par22+par33;
   KSGC(i)=par11;
end
%%计算pc值，找到初始群C，吸收节点，循环得到全部的群落Cs

mix=adjMat;
N=size(mix,1);  %节点数
degree=sum(mix,2)'; %按行求和，得节点的度
avgdeg=sum(degree)/N;  %平均度
for i=1:N
    B1=find(mix(i,:));%和第i个节点相连的节点索引
    len1=length(B1);%len1表示i的邻居数
    ddin(i)=0;
    for j=1:len1
        B2=find(mix(B1(j),:));%i的邻居之一j的邻居
        jiaoji=intersect(B1,B2);%i与j的邻居交集
        ddin(i)=length(jiaoji)+ddin(i);
    end
    pc(i)=degree(i)*sum(ddin(i));
end
core=pc;
C={};
ii=1;
totlenc=0;
num=0;
threshold=0;
while threshold<=0.6  %达到网络节点总数的60%停止

    
    %[a,b]=sort(core,'descend'); %a-降序序列 b-节点索引
    [pcmax,index]=max(core);
    biggest=find(core==pcmax);%找到pcmax的节点
    aaa=biggest(find(degree(biggest)==max(degree(biggest))));%找到pc值最大的里面度最大的
    len0=length(aaa);
    initial_set(ii)=aaa(ceil(rand(1)*len0));%随机取得最大pc值中的一个作为起始节点%%%%
    num=num+1;
    if num/N<=0.6
        %initial_set(ii)=aaa(1);
        C{ii}=initial_set(ii);
        neighbor_of_i=find(mix(initial_set(ii),:))%初始节点的邻居都找出来
        zz=find(pc(neighbor_of_i)>=pc(initial_set(ii))/2);%找到pc值大于1/2pcmax的节点
        pc_of_zz=pc(neighbor_of_i(zz));
        togther0=[pc_of_zz;neighbor_of_i(zz)]';
        togtherr0=sortrows(togther0)';
        for j=length(zz):-1:1
            C{ii}=[C{ii},togtherr0(2,j)];
            num=num+1;
            if num/N>0.6
                C{ii} = C{ii}(1:end-1);
                break
            end
        end

        % C{ii}=neighbor_of_i(zz);%找到节点i的内部链接的各个节点,组成子图C
        % num=num+length(C{ii});
        %C{ii}=[C{ii},initial_set(ii)];
        if num/N<=0.6
            first_neighbor_of_C=neighbor_of_node_set(mix,C{ii});
            degree_of_first_neighbor_of_C=degree(first_neighbor_of_C);
            togther=[degree_of_first_neighbor_of_C;first_neighbor_of_C]';
            togtherr=sortrows(togther)';
            
            %for j=length(first_neighbor_of_C):-1:1
            for j=1:length(first_neighbor_of_C)
                ddd=neighbor(mix,togtherr(2,j));
                jiaoji=intersect(ddd,C{ii});
                chaji=setdiff(ddd,jiaoji);
                if length(jiaoji)>length(chaji)
                    C{ii}=[C{ii},togtherr(2,j)];
                    num=num+1;
                    if num/N>0.6
                        C{ii} = C{ii}(1:end-1);
                        break
                    end
                end
            end
        else
            break
        end
        if num/N<=0.6
            second_neighbor_of_C=neighbor_of_node_set(mix,C{ii});
            degree_of_second_neighbor_of_C=degree(second_neighbor_of_C);
            togther=[degree_of_second_neighbor_of_C;second_neighbor_of_C]';
            togtherr=sortrows(togther)';
            %for j=length(second_neighbor_of_C):-1:1
            for j=1:length(second_neighbor_of_C)
                ddd=neighbor(mix,togtherr(2,j));
                jiaoji=intersect(ddd,C{ii});
                chaji=setdiff(ddd,jiaoji);
                if length(jiaoji)>length(chaji)
                    C{ii}=[C{ii},togtherr(2,j)];
                    num=num+1;
                    if num/N>0.6
                        C{ii} = C{ii}(1:end-1);
                        break
                    end
                end
        
            end
        else
            break
        end
        if num/N<=0.6
            third_neighbor_of_C=neighbor_of_node_set(mix,C{ii});
            degree_of_third_neighbor_of_C=degree(third_neighbor_of_C);
            togther=[degree_of_third_neighbor_of_C;third_neighbor_of_C]';
            togtherr=sortrows(togther)';
            %for j=length(third_neighbor_of_C):-1:1
            for j=1:length(third_neighbor_of_C)
                ddd=neighbor(mix,togtherr(2,j));
                jiaoji=intersect(ddd,C{ii});
                chaji=setdiff(ddd,jiaoji);
                if length(jiaoji)>length(chaji)
                    C{ii}=[C{ii},togtherr(2,j)];
                    num=num+1;
                    if num/N>0.6
                        C{ii} = C{ii}(1:end-1);
                        break
                    end
                end
        
            end
        else
            break
        end
    
        %C_number_of_initial_i(i)=length(C{i});
        linshi3=C{ii};
        lenc=length(linshi3);
        totlenc=totlenc+lenc;
        linshi4=linshi3(find(degree(linshi3)==max(degree(linshi3))));
        initial_set(ii)=linshi4(ceil(rand(1)*length(linshi4)));
        %initial_set(ii)=linshi4(1);
    
        mix(C{ii},:)=0;
        mix(:,C{ii})=0;
        degree=sum(mix);
        MN=sparse(mix);
        for i=1:N
            B1=find(mix(i,:));%和第i个节点相连的节点索引
            len1=length(B1);%len1表示i的邻居数
            ddin(i)=0;
            for j=1:len1
                B2=find(mix(B1(j),:));%i的邻居之一j的邻居
                jiaoji=intersect(B1,B2);%i与j的邻居交集
                ddin(i)=length(jiaoji)+ddin(i);
            end
            core(i)=degree(i)*sum(ddin(i));
        end
        
        
        ii=ii+1;
        threshold=totlenc/N;    
    else
        initial_set=initial_set(1:end-1);
        break
    end
end



%%节点局部影响计算
kcore=core_numbers(sparse(adjMat))';
LI=[];
for i=1:N
    B3=find(adjMat(i,:));%和第i个节点相连的节点索引
    len2=length(B3);%len3表示i的邻居数
    ks=0;
    for j=1:len2
        ks=ks+kcore(B3(j));
    end
    LI(end+1)=ks;
end
shortest_distances = dijkstra_for_targets(adjMat,initial_set);%找到全局关键节点到其他节点的最短距离

%%

sd=shortest_distances';
for i=1:N
    GAC(i)=0;
    for j=1:length(initial_set)
        GAC(i)=LI(i)*LI(initial_set(j))/(2^sd(i,j))+GAC(i);
    end
end

average_i=sir_main(adjMat,1000,0.054) ;%传播概率 每个节点的平均感染数0.047

save('email.mat', 'd', 'h','p','MDD_core','kcore','neighbor_Core','KSGC','GAC','DC_stc','average_i');

average_iall={};
w = waitbar(0, 'Please wait...'); % 初始化进度条

for i = 1:15
    waitbar(i/15,w , sprintf('Progress: %d %%', floor(i/15*100))); % 更新进度条
    disp(i);
    average_iall{i} = sir_main(adjMat, 1000, 0.01*i);
end

close(w); % 运行结束后关闭进度条

%a00_hefei1=[];a00_h=[];a00_pagerank=[];a00_M_core=[];a00_neighbor_core=[];a00_degree=[];a0_core=[];a00_KSGC=[];a00_h=[];
for i=1:15
    a00_degree(i)=corr(average_iall{i}',d','type','Kendall');%度中心性
    a00_h_index(i)=corr(average_iall{i}',h,'type','Kendall');%H指数
    a00_pagerank(i)=corr(average_iall{i}',p,'type','Kendall');%PageRank中心性
    a00_MDD_core(i)=corr(average_iall{i}',MDD_core','type','Kendall');%混合度分解法MDD
    a00_GAC(i)=corr(average_iall{i}',GAC','type','Kendall');
    a00_kcore(i)=corr(average_iall{i}',kcore','type','Kendall');%k壳分解法
    a00_neighbor_core(i)=corr(average_iall{i}',neighbor_Core','type','Kendall');%邻域核心度中心性（NCC）
    a00_KSGC(i)=corr(average_iall{i}',KSGC','type','Kendall');%基于 K 核分解法的引力模型改进方法KSGC
    a00_TLGM(i)=corr(average_iall{i}',DC_stc,'type','Kendall');
end

save('email_sir1_15_Kendall.mat', 'average_iall','a00_degree', 'a00_h_index','a00_pagerank','a00_MDD_core','a00_GAC','a00_kcore','a00_neighbor_core','a00_KSGC','a00_TLGM');
%is_symmetric = issymmetric(all_resitace_distance)
M=zeros(1, 8);
M(1)=calculate_monotonicity(d);
M(2) = calculate_monotonicity(p);
M(3)=calculate_monotonicity(h);
M(4) = calculate_monotonicity(kcore);
M(5)=calculate_monotonicity(MDD_core);
M(6) = calculate_monotonicity(neighbor_Core);
M(7)=calculate_monotonicity(KSGC);
M(8) = calculate_monotonicity(DC_stc);
ccdf=cell(8, 1);
ccdf{1}=CCDF(d);
ccdf{2}=CCDF(p);
ccdf{3}=CCDF(h);
ccdf{4}=CCDF(kcore);
ccdf{5}=CCDF(MDD_core);
ccdf{6}=CCDF(neighbor_Core);
ccdf{7}=CCDF(KSGC);
ccdf{8}=CCDF(DC_stc);

x = 0.01:0.01:0.15;  % x轴数据
% 调整整个图的尺寸
fig = gcf;
fig.Position = [100, 100, 1000, 600]; % [left, bottom, width, height]
% 为每种算法绘制曲线，并设置颜色和标记
plot(x, a00_GAC, 's-', 'Color', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0], 'DisplayName', 'GLC', 'MarkerSize', 6); % 黑色
hold on;
plot(x, a00_h_index, 'o-', 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'H', 'MarkerSize', 6); % 红色
plot(x, a00_pagerank, '^-', 'Color', [0, 0, 1], 'MarkerFaceColor', [0, 0, 1], 'DisplayName', 'PRC', 'MarkerSize', 6); % 蓝色
plot(x, a00_MDD_core, 'v-', 'Color', [1, 0, 1], 'MarkerFaceColor', [1, 0, 1], 'DisplayName', 'MDD', 'MarkerSize', 6); % 暗洋红色
plot(x, a00_neighbor_core, 'd-', 'Color', [0, 128/255, 0], 'MarkerFaceColor',[0, 128/255, 0], 'DisplayName', 'NCC', 'MarkerSize', 6); % 暗黄色
plot(x, a00_degree, '<-', 'Color', [0, 0, 128/255], 'MarkerFaceColor', [0, 0, 128/255], 'DisplayName', 'DC', 'MarkerSize', 6); % 暗红色
plot(x, a00_kcore, '>-', 'Color', [128/255, 0, 1], 'MarkerFaceColor', [128/255, 0, 1], 'DisplayName', 'KS', 'MarkerSize', 6); % 深灰色
plot(x, a00_KSGC, 'h-', 'Color', [128/255, 0, 128/255], 'MarkerFaceColor', [128/255, 0, 128/255], 'DisplayName', 'KSGC', 'MarkerSize', 6); % 深粉红色
% 添加 a00_DC_stc 曲线（青绿色）
plot(x, a00_TLGM, 'p-', 'Color', [0, 1, 1], 'MarkerFaceColor', [0, 1, 1], 'DisplayName', 'DC_stc', 'MarkerSize', 6); % 青色
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