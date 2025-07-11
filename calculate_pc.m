function pc = calculate_pc(mixedsig, core, EKC_values)
    % 获取节点总数
    N = length(core);
    % 初始化 pc 矩阵
    pc = zeros(N, 1);
    % 遍历每个节点
    for i = 1:N
        % 获取节点 i 的一阶邻居
        neighbors = find(mixedsig(i, :) == 1);
        % 计算一阶邻居的 EKC 值之和
        sum_EKC_neighbors = sum(EKC_values(neighbors));
        % 计算节点 i 的 pc 值
        pc(i) = core(i) * sum_EKC_neighbors;
    end
end
    