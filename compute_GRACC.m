function GRACC = compute_GRACC(adjMat)
    % 计算初始网络的聚类系数
    CC_values = clustering_coefficient(adjMat);  % 计算每个节点的聚类系数
    avg_CC = mean(CC_values);  % 初始网络的平均聚类系数
    
    N = size(adjMat, 1);
    GRACC = zeros(1, N);
    
    for v = 1:N
        % 移除节点 v 后重新计算网络的聚类系数
        adjMat_v = adjMat;
        adjMat_v(v, :) = 0;  % 移除节点 v 的所有出边
        adjMat_v(:, v) = 0;  % 移除节点 v 的所有入边
        
        % 重新计算移除节点 v 后的聚类系数
        CC_values_v = clustering_coefficient(adjMat_v);
        
        % 在计算移除节点后的平均聚类系数时，除以剩余节点数 N - 1
        avg_CC_v = sum(CC_values_v) / (N - 1);
        
        % 计算 GRACC
        GRACC(v) = abs(avg_CC_v - avg_CC) / avg_CC;
    end
end