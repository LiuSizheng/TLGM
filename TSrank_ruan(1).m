clear;
    % 加载数据集并处理成邻接矩阵
    A = load('jazz.txt');
    % A = A + 1;
    r = 0.2;
    TT = A(:, 1:2);
    mixedsig = zeros(max(max(TT)));
    N = size(mixedsig, 1);
    len = length(TT);
    for i = 1:len
        mixedsig(TT(i, 1), TT(i, 2)) = 1;
        mixedsig(TT(i, 2), TT(i, 1)) = 1;
    end
    kk = sparse(mixedsig);
    [a, b] = components(kk);
    [B] = largestcomponent(mixedsig);
    mixedsig = mixedsig(B, B);

    % 计算结构洞约束系数
    cc = burt(mixedsig);
    % 计算 MDD 得到 Ksi 值
    M_core = MDD(mixedsig);

    % 计算节点数量
    numNodes = size(mixedsig, 1);
    % 初始化存储结果的向量
    TSM_values = zeros(numNodes, 1);

    % 计算度中心性
    DC = sum(mixedsig, 2);
    % 计算邻接度
    Q = zeros(numNodes, 1);
    for i = 1:numNodes
        neighbors = find(mixedsig(i, :));
        for j = 1:length(neighbors)
            Q(i) = Q(i) + DC(neighbors(j));
        end
    end

 

    % 计算 K 壳值
    KS = M_core;

   

    % 计算节点信息熵
    I = DC / sum(DC);
    e = zeros(numNodes, 1);
  

  % 计算 Tsallis 熵相关参数
    IC = zeros(numNodes, 1); % 存储每个节点的 IC 值
    for i = 1:numNodes
        % 找到节点 i 的邻居
        neighbors = find(mixedsig(i, :)); 
        Pij = zeros(length(neighbors), 1);
        for j = 1:length(neighbors)
            % 根据结构洞约束系数计算 Pij
            Pij(j) = (1 - cc(neighbors(j))) / sum(1 - cc(neighbors)); 
        end
        % Ksi 为节点 i 的 K 壳值与所有节点 K 壳值之和的比值
        Ksi = KS(i) / sum(KS); 
        qi = 1 + Ksi; 
        % 计算 Tsallis 熵 T
        T = (1 - sum(Pij.^qi)) / (qi - 1); 
        % 结合结构洞约束系数计算 IC
        IC(i) = (1 - cc(i)) * T; 
    end
    
    % 计算 Cnc 和 TSM 值
    for i = 1:numNodes
        % 找到节点 i 的邻居
        neighbors = find(mixedsig(i, :)); 
        Cnc = 0;
        for j = 1:length(neighbors)
            % 计算邻域核心
            Cnc = Cnc + IC(neighbors(j)); 
        end
        % 计算 TSM 值
        TSM_values(i) = Cnc; 
    end

    % 输出最终的 TSM 节点重要性评估值
    disp(TSM_values); 