function cc = closeness(adj_matrix)
    N = size(adj_matrix, 1); % 节点数量
    closeness = zeros(N, 1); % 初始化接近中心性向量
    all_nodes = 1:N;
    
    for v = 1:N
        % 使用 Dijkstra 算法计算节点 v 到所有其他节点的最短路径距离
        distances = dijkstra_for_targets(adj_matrix, v);
        
        % 计算可达节点的数量（包括自身）
        reachable = sum(distances < inf);
        
        if reachable > 1 % 至少能到达其他节点
            % 总的距离和（排除自身）
            total_distance = sum(distances(distances < inf)) - distances(v);
            % 计算接近中心性
            cc(v) = (reachable - 1) / total_distance;
        else
            cc(v) = 0; % 不可达的情况，接近中心性定义为 0
        end
    end
end