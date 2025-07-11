function algebraic_distance = calculate_algebraic_distance1(mixedsig)
    % 计算度矩阵
    D = diag(sum(mixedsig));
    % 计算拉普拉斯矩阵
    L = D - mixedsig;
    % 对拉普拉斯矩阵进行特征分解
    [U, Lambda] = eig(L);
    % 获取节点数量
    n = size(mixedsig, 1);
    % 初始化代数距离矩阵
    algebraic_distance = zeros(n, n);

    % 计算代数距离
    for i = 1:n
        % 找出 1 跳邻居
        one_hop_neighbors = find(mixedsig(i, :));
        % 找出 2 跳邻居
        two_hop_neighbors = [];
        for neighbor = one_hop_neighbors
            two_hop = find(mixedsig(neighbor, :));
            two_hop_neighbors = [two_hop_neighbors, two_hop];
        end
        two_hop_neighbors = unique(two_hop_neighbors);
        % 去除自身和 1 跳邻居
        two_hop_neighbors = setdiff(two_hop_neighbors, [one_hop_neighbors, i]);
        all_neighbors = unique([one_hop_neighbors, two_hop_neighbors]);

        for j = all_neighbors
            for k = 1:n
                if Lambda(k, k) ~= 0
                    algebraic_distance(i, j) = algebraic_distance(i, j) + (U(i, k) - U(j, k))^2 / Lambda(k, k);
                end
            end
            algebraic_distance(i, j) = sqrt(algebraic_distance(i, j));
        end
    end
end

