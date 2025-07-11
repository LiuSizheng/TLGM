%计算每一个节点的EWD
function ewd=EWD(mixedsig) 


mix=mixedsig;


%% RMD 分解算法实现
N = size(mix, 1); % 当前图的节点数
CurrentIter = 1; % 当前迭代次数
Iter = zeros(N, 1); % 记录每个节点的RMD索引

% 创建一个用于标记哪些节点已经被删除的数组
deletedNodes = false(N, 1);

while N > 0
    % 计算每个节点的度数（忽略已经被删除的节点）
    degree = sum(mix, 2);
    degree(deletedNodes) = inf; % 对已经删除的节点度数设为无穷大，避免它们被再次选择

    % 找到最小度数
    md = min(degree(degree < inf)); % 找到度数中最小的值，但忽略无穷大部分

    % 存储将要删除的节点
    DelArray = find(degree == md);

    % 给这些节点分配当前的迭代索引
    Iter(DelArray) = CurrentIter;

    % 从邻接矩阵中删除这些节点
    for v = DelArray'
        mix(v, :) = 0;
        mix(:, v) = 0;
        deletedNodes(v) = true; % 标记该节点为已删除
    end

    % 更新节点数和当前迭代次数
    N = N - length(DelArray);
    CurrentIter = CurrentIter + 1;
end

% 计算最大RMD索引
MaxIter = max(Iter);
%计算WD
n = size(mixedsig, 1); % 当前图的节点数
WD = zeros(n, 1); 

for i=1:n
    for j=1:n
        if mixedsig(i,j)==1;
            WD(i)=WD(i)+Iter(j)/MaxIter;
        end 
    end
end
%计算EWD
ewd=zeros(n, 1); %
for k=1:n
    for p=1:n
        if mixedsig(k,p)==1;
            ewd(k)=ewd(k)+WD(p);
        end 
    end
end





