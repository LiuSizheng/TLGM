function distances = compute_resistance_distance(A, i)
    % 输入：
    % A - 邻接矩阵 (无向无权图，对称矩阵)
    % i - 目标节点编号
    % 输出：
    % distances - 3×m 的二维数组，m为三阶邻域内节点数
    %           第一行：与一阶邻居的有效距离
    %           第二行：与二阶邻居的有效距离
    %           第三行：与三阶邻居的有效距离
    
    n = size(A, 1);  % 图的节点数
    
    % 第一步：找到一阶邻居
    first_order = find(A(i,:) == 1);
    
    % 第二步：找到二阶邻居
    second_order = [];
    for j = first_order
        neighbors = find(A(j,:) == 1);
        second_order = union(second_order, neighbors);
    end
    second_order = setdiff(second_order, [i first_order]);
    
    % 第三步：找到三阶邻居
    third_order = [];
    for j = second_order
        neighbors = find(A(j,:) == 1);
        third_order = union(third_order, neighbors);
    end
    third_order = setdiff(third_order, [i first_order second_order]);
    
    % 所有三阶以内邻居（包括一阶、二阶、三阶）
    all_neighbors = union(union(first_order, second_order), third_order);
    
    % 初始化距离数组：3行，列数为所有邻居数量
    m = length(all_neighbors);
    distances = zeros(3, m);  % 3行分别对应1阶、2阶、3阶邻居距离
    
    % 创建邻居到索引的映射
    neighbor_indices = containers.Map(all_neighbors, 1:m);
    
    % 对每个邻居计算有效距离
    for j = all_neighbors
        idx = neighbor_indices(j);  % 当前邻居在distances中的列索引
        resist_sum = 0;
        
        % 一阶路径 (直接连接)
        if A(i,j) == 1
            resist_sum = resist_sum + 1/1;  % 1步路径电阻倒数
        end
        
        % 二阶路径
        common_neighbors_2 = intersect(first_order, find(A(j,:) == 1));
        num_paths_2 = length(common_neighbors_2);
        resist_sum = resist_sum + num_paths_2 * (1/2);  % 2步路径电阻倒数
        
        % 三阶路径
        num_paths_3 = 0;
        for k = second_order
            if A(k,j) == 1  % k是j的一阶邻居
                common_neighbors_ik = intersect(first_order, find(A(k,:) == 1));
                num_paths_ik = length(common_neighbors_ik);
                num_paths_3 = num_paths_3 + num_paths_ik;
            end
        end
        resist_sum = resist_sum + num_paths_3 * (1/3);  % 3步路径电阻倒数
        
        % 根据邻居的阶数存储有效距离
        if resist_sum > 0
            if any(first_order == j)
                distances(1, idx) = 1/resist_sum;  % 一阶邻居
            elseif any(second_order == j)
                distances(2, idx) = 1/resist_sum;  % 二阶邻居
            elseif any(third_order == j)
                distances(3, idx) = 1/resist_sum;  % 三阶邻居
            end
        end
    end
end

