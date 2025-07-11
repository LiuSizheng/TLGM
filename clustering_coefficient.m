function CLC = clustering_coefficient(adjMat)
    % 计算每个节点的聚类系数
    N = size(adjMat, 1);
    CLC = zeros(1, N);
    
    for i = 1:N
        neighbors = find(adjMat(i, :) == 1);
        k_i = length(neighbors);
        
        if k_i > 1
            subgraph = adjMat(neighbors, neighbors);
            links = sum(subgraph(:)) / 2;
            CLC(i) = 2 * links / (k_i * (k_i - 1));
        else
            CLC(i) = 0;
        end
    end
end