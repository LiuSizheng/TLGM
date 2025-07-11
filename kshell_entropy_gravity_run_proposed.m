function G_values = kshell_entropy_gravity_run_proposed(mixedsig, core, N,result)

% 存储每个节点的 EKC 值
EKC_values = zeros(N, 1);
for i = 1:N
    % 计算一阶邻居的 k-壳熵值 EK1
    first_neighbors = find(mixedsig(i, :));
    if ~isempty(first_neighbors)
        k_shell_values_1st = core(first_neighbors);
        K_sum1 = sum(k_shell_values_1st);
        EK1 = -sum((k_shell_values_1st / K_sum1).* log(k_shell_values_1st / K_sum1));
        % 计算二阶邻居的 k-壳熵值 EK2
        second_neighbors = [];
        for j = first_neighbors
            second_neighbors = [second_neighbors, find(mixedsig(j, :))];
        end
        second_neighbors = unique(second_neighbors);
        % 排除节点 i 本身和一阶邻居
        second_neighbors = setdiff(second_neighbors, [i, first_neighbors]);
        k_shell_values_2nd = core(second_neighbors);
        K_sum2 = sum(k_shell_values_2nd);
        EK2 = -sum((k_shell_values_2nd / K_sum2).* log(k_shell_values_2nd / K_sum2));
        % 计算 lambda
        lambda = 1 / (1 + exp(-(EK2 - EK1)));
        % 计算节点的 k-壳多样性 EKC
        EKC_values(i) = EK1 + lambda * EK2;
    else
        EKC_values(i) = 0; % 对于孤立节点，EKC 值设为 0
    end
end
KEKC = calculate_pc(mixedsig, core, EKC_values);


% 计算 G(i)，将一阶、二阶、三阶邻居分开计算
G_values = zeros(N, 1);
for i = 1:N
    % 一阶邻居
    first_neighbors = find(mixedsig(i, :));
    G_first = 0;
    for j = first_neighbors
        % 使用谱距离替换原来的距离
        dist = result(i, j);
        G_first = G_first + (KEKC(i) * KEKC(j)) / ((1+dist)^1);
    end
    
    % 二阶邻居
    second_neighbors = [];
    for j = first_neighbors
        second_neighbors = [second_neighbors, find(mixedsig(j, :))];
    end
    second_neighbors = unique(second_neighbors);
    second_neighbors = setdiff(second_neighbors, [i, first_neighbors]); % 去除自身和一阶邻居
    G_second = 0;
    for j = second_neighbors
        % 使用谱距离替换原来的距离
        dist = result(i, j);
        G_second = G_second + (KEKC(i) * KEKC(j)) / ((1+dist)^2);
    end
    
%     % 三阶邻居
%     third_neighbors = [];
%     for j = second_neighbors
%         third_neighbors = [third_neighbors, find(A(j, :))];
%     end
%     third_neighbors = unique(third_neighbors);
%     third_neighbors = setdiff(third_neighbors, [i, first_neighbors, second_neighbors]); % 去除自身、一阶和二阶邻居
%     G_third = 0;
%     for j = third_neighbors
%         % 使用谱距离替换原来的距离
%         dist = result(i, j);
%         G_third = G_third + (KEKC(i) * KEKC(j)) / ((1+dist)^3);
%     end
    
    G_values(i) = G_first + G_second+0;
end
end