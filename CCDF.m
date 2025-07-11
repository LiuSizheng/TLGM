function ccdf=CCDF(data)
    % PLOT_CCDF 计算并绘制数据的互补累积分布函数 (CCDF)
    % 输入:
    %   data - 数据向量（例如节点排名）
    % 输出:
    %   无直接输出，生成 CCDF 图
    
    % 确保数据是列向量
    data = data(:);
    
    % 去重并排序数据
    unique_values = sort(unique(data), 'descend');
    n = length(data); % 数据总数
    
    % 计算每个唯一值的 CCDF
    ccdf = zeros(size(data));
    for i = 1:length(unique_values)
        ccdf(i) = (n-sum(data >= unique_values(i))) / n;
    end
  
    
end