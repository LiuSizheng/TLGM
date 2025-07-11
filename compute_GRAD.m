function GRAD = compute_GRAD(adjMat)
    % 计算网络的平均度
    degree = sum(adjMat, 2);  % 计算每个节点的度
    avg_degree = mean(degree);  % 计算平均度
    
    N = size(adjMat, 1);
    GRAD = zeros(1, N);
    
    for v = 1:N
        % 移除节点 v 并重新计算平均度
        adjMat_v = adjMat;
        adjMat_v(v, :) = 0;
        adjMat_v(:, v) = 0;
        
        % 计算移除节点后的平均度
        degree_v = sum(adjMat_v, 2);
        avg_degree_v = sum(degree_v)/(N-1);
        
        % 计算 GRAD
        GRAD(v) = abs(avg_degree_v - avg_degree) / avg_degree;
    end
end