function neighbors = find_local_neighbors(adjMat, v, L)
    % 寻找节点 v 的 L 阶邻居
    N = size(adjMat, 1);
    neighbors = false(1, N);
    neighbors(v) = true;
    
    for i = 1:L
        neighbors = neighbors | any(adjMat(neighbors, :), 1);
    end
    
    neighbors = find(neighbors);
end