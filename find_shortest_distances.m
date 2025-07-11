function shortest_distances = find_shortest_distances(N, M, adj_matrix)
    % 初始化距离矩阵为无穷大
    distance_matrix = inf(N, N);
    
    % 对于邻接矩阵中的每一个1，设置对应的距离为1
    for i = 1:N
        for j = 1:N
            if adj_matrix(i, j) == 1
                distance_matrix(i, j) = 1;
            end
        end
    end
    
    % 设置对角线元素为0
    for i = 1:N
        distance_matrix(i, i) = 0;
    end
    
    % Floyd-Warshall算法
    for k = 1:N
        for i = 1:N
            for j = 1:N
                distance_matrix(i, j) = min(distance_matrix(i, j), distance_matrix(i, k) + distance_matrix(k, j));
            end
        end
    end
    
    shortest_distances = distance_matrix;
end