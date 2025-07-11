% 加载原始 .mat 文件内容到结构体
data = load('netscience_sir1_26_Kendall.mat');

% 删除指定的变量
if isfield(data, 'a00_EKC_values')
    data = rmfield(data, 'a00_EKC_values');
end
if isfield(data, 'a00_neighbor_EKC_values')
    data = rmfield(data, 'a00_neighbor_EKC_values');
end
if isfield(data, 'a00_extend_neighbor_EKC')
    data = rmfield(data, 'a00_extend_neighbor_EKC');
end

% 保存更新后的数据回原文件
save('netscience_sir1_26_Kendall.mat', '-struct', 'data');