function LRACC = compute_LRACC(adjMat, L)
    N = size(adjMat, 1);
    LRACC = zeros(1, N);
    
    % 计算每个节点的聚类系数
    
    
    for v = 1:N
        % 获取节点 v 的局部邻居（距离不超过 L 的节点）
        neighbors = find_local_neighbors(adjMat, v, L);
        
        % 创建局部邻域的子图
        subgraph = adjMat(neighbors, neighbors);
        CC_values = clustering_coefficient(subgraph);
        % 计算移除节点后的局部聚类系数
        subgraph_v = subgraph;
        v_index = find(neighbors == v);
        subgraph_v(v_index, :) = 0;
        subgraph_v(:, v_index) = 0;
        
        % 计算聚类系数的平均变化
        CC_values_v = clustering_coefficient(subgraph_v);
        avg_CC = mean(CC_values);
        avg_CC_v = sum(CC_values_v) / (numel(neighbors) - 1);  % 除以剩余的邻居节点数量
        
        % 计算 LRACC
        LRACC(v) = abs(avg_CC_v - avg_CC) / avg_CC;
    end
end
