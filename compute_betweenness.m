function betweenness = compute_betweenness(adj_matrix)
    % 计算无向无权图的节点介数中心性
    N = size(adj_matrix, 1); % 获取节点数量
    betweenness = zeros(N, 1); % 初始化所有节点的介数中心性为零

    for s = 1:N
        % 初始化变量
        stack = []; % 用于回溯的栈
        P = cell(N, 1); % 前驱节点列表
        sigma = zeros(N, 1); % 最短路径数量
        sigma(s) = 1; % 源节点的最短路径数量初始化为 1
        d = -1 * ones(N, 1); % 距离数组，初始为 -1 表示未访问
        d(s) = 0; % 源节点到自身的距离为 0
        Q = [s]; % 初始化队列，将源节点加入队列

        % 广度优先搜索 (BFS)
        while ~isempty(Q)
            v = Q(1); % 取队首元素
            Q(1) = []; % 将队首出列
            stack = [stack, v]; % 将节点压入栈中，以便回溯

            neighbors = find(adj_matrix(v, :) > 0); % 找到邻居节点
            for w = neighbors
                % 路径发现
                if d(w) < 0 % 节点 w 第一次被访问
                    d(w) = d(v) + 1; % 更新距离
                    Q = [Q, w]; % 将 w 加入队列
                end
                % 计算最短路径数量
                if d(w) == d(v) + 1 % 如果 w 在最短路径上
                    sigma(w) = sigma(w) + sigma(v); % 更新 w 的路径数量
                    P{w} = [P{w}, v]; % 将 v 记录为 w 的前驱节点
                end
            end
        end

        % 回溯计算依赖度 delta
        delta = zeros(N, 1); % 初始化依赖度
        while ~isempty(stack)
            w = stack(end); % 取栈顶元素
            stack(end) = []; % 出栈

            for v = P{w} % 遍历 w 的所有前驱节点
                delta(v) = delta(v) + (sigma(v) / sigma(w)) * (1 + delta(w)); % 更新依赖度
            end

            if w ~= s % 排除源节点自身
                betweenness(w) = betweenness(w) + delta(w); % 累加到介数中心性
            end
        end
    end

    % 无向图的中心性值需要除以 2
    betweenness = betweenness / 2;
end
