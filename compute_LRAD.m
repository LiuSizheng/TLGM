function LRAD = compute_LRAD(adjMat, L)
    N = size(adjMat, 1);
    LRAD = zeros(1, N);
    
    for v = 1:N
        % 获取节点 v 的局部邻居（距离不超过 L 的节点）
        neighbors = find_local_neighbors(adjMat, v, L);
        
        % 创建局部邻域的子图
        subgraph = adjMat(neighbors, neighbors);
        
        % 计算移除节点前局部邻域的平均度
        degree_values = sum(subgraph, 2);  % 计算子图中的度
        avg_degree = mean(degree_values);  % 移除前的平均度
        
        % 移除节点 v 的影响
        subgraph_v = subgraph;
        v_index = find(neighbors == v);
        subgraph_v(v_index, :) = 0;
        subgraph_v(:, v_index) = 0;
        
        % 重新计算移除节点 v 后的平均度
        degree_values_v = sum(subgraph_v, 2);
        avg_degree_v = sum(degree_values_v) / (numel(neighbors) - 1);  % 除以剩余的邻居节点数
        
        % 计算 LRAD
        LRAD(v) = abs(avg_degree_v - avg_degree) / avg_degree;
    end
end
