function kendalls_tau = kendalls_tau(x, y)
    % 检查输入数组长度是否一致
    if length(x) ~= length(y)
        error('两个数组的长度必须相等');
    end

    % 初始化一致对和不一致对计数
    n_cp = 0;  % 一致对数量
    n_dp = 0;  % 不一致对数量

    % 遍历所有可能的对 (i, j) 且 i < j
    n = length(x);
    for i = 1:n-1
        for j = i+1:n
            % 计算 (x_i - x_j) * (y_i - y_j)
            sign_product = (x(i) - x(j)) * (y(i) - y(j));
            if sign_product > 0
                n_cp = n_cp + 1;  % 一致对
            elseif sign_product < 0
                n_dp = n_dp + 1;  % 不一致对
            end
        end
    end

    % 总对数 tp = n * (n - 1) / 2
    tp = n * (n - 1) / 2;

    % 计算 Kendall's Tau 相关系数
    kendalls_tau = (n_cp - n_dp) / tp;
end
