function ISBC = compute_ISBC(adj_matrix)
    % 获取节点数量
    N = size(adj_matrix, 1);
    
    % 计算每个节点的介数中心性 (BC)
    BC = compute_betweenness(adj_matrix)*2/((N-1)*(N-2)) ;
    
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
    
    % 使用dijkstra_for_targets计算每个节点的最短路径距离
    target_nodes = 1:N;
    shortest_distances = dijkstra_for_targets(adj_matrix, target_nodes);
    
    % 计算每个节点的 ISBC 值
    ISBC = zeros(1, N); % 初始化 ISBC 值
    for v = 1:N
        sum_val = 0;
        for u = 1:N
            if u ~= v
                d_uv = shortest_distances(v, u); % 从最短路径距离矩阵中获取节点 v 和 u 的距离
                if d_uv > 0 % 避免除以零
                    sum_val = sum_val + sqrt(BC(u)) / d_uv; %sum_val = sum_val + sqrt(BC(u)) / d_uv;
                end
            end
        end
        ISBC(v) = (ISC(v) / N) * sum_val;
    end
end
