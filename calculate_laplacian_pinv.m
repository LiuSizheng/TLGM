function L_pinv = calculate_laplacian_pinv(mixedsig)
    % 计算度矩阵
    D = diag(sum(mixedsig));
    % 计算拉普拉斯矩阵
    L = D - mixedsig;
    % 计算拉普拉斯矩阵的伪逆矩阵（摩尔-彭罗斯伪逆）
    L_pinv = pinv(L);
end