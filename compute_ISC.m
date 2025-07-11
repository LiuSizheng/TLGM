function ISC = compute_ISC(adj_matrix)
    % 获取节点数量
    N = size(adj_matrix, 1);
    
    % 计算每个节点的隔离中心性 (ISC)
    ISC = zeros(1, N); % 初始化隔离中心性为零
    degrees = sum(adj_matrix, 2); % 计算每个节点的度数
    min_degree = min(sum(adj_matrix, 2)); % 计算网络中的最小度数

    
    for v = 1:N
        neighbors = find(adj_matrix(v, :) > 0); % 获取邻居节点
        degree_v = degrees(v); % 计算节点 v 的度
        if isempty(neighbors) % 如果没有邻居，则 ISC 为 0
            ISC(v) = 0;
            continue;
        end
       
        neighbor_degrees = degrees(neighbors);
      
        % 计算 ISC: 统计邻居中度数为最小的节点个数
        num_min_degree_neighbors = sum(neighbor_degrees == min_degree);
        ISC(v) = num_min_degree_neighbors * degree_v; % ISC 计算
    end
end