function diameter = calculate_diameter(adjMat) %计算网络直径
    N = size(adjMat, 1);
    max_distance = 0;
    
    for i = 1:N
        distances = dijkstra_for_targets(adjMat, i); % 使用dijkstra_for_targets计算i节点到其他节点的最短路径
        max_distance = max(max_distance, max(distances(distances < inf))); % 更新最大距离，忽略无穷大值
    end
    
    diameter = max_distance;
end
