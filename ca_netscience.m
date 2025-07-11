%%
%%读取网络文件，找到最大连通分量
clear;
net=load('ca-netscience.txt'); %读取网络文件test0425
adjMat=zeros(max(max(net))); %创建网络节点数量大小的邻接矩阵，初始为0
%N=size(adjMat,1); %获取行数
len=length(net); %获取net行数
for i =1:len   %邻接赋值为1
    adjMat(net(i,1),net(i,2))=1;
    adjMat(net(i,2),net(i,1))=1;
end
%spMat=sparse(adjMat);%稀疏矩阵
%[a b]=components(spMat);%找到连通分量 a-连通分量序号 b-分量节点数
[B]=largestcomponent(adjMat); %B返回最大连通分量包含的节点序号
adjMat=adjMat(B,B); %得到最大连通分量的邻接矩阵
mix=adjMat;
N=size(mix,1);  %节点数
degree=sum(mix,2)'; %按行求和，得节点的度
avgdeg=sum(degree)/N;  %平均度
DC = degreea( adjMat ); %度中心性
BC = compute_betweenness(adjMat)'*2/((N-1)*(N-2));
CC=closeness(adjMat);
ISBC = compute_ISBC(adjMat);
ISC=compute_ISC(adjMat);

ISBC_DC_K=kendalls_tau(ISBC,DC);
ISBC_BC_K=kendalls_tau(ISBC,BC);
ISBC_CC_K=kendalls_tau(ISBC,CC);
ISBC_ISC_K=kendalls_tau(ISBC,ISC);

% ISBC_DC_K=corr(ISBC',DC','type','Kendall');%度中心性
% ISBC_BC_K=corr(ISBC',BC','type','Kendall');%H指数
% ISBC_CC_K=corr(ISBC',CC','type','Kendall');%PageRank中心性
% ISBC_ISC_K=corr(ISBC',ISC','type','Kendall');%混合度分解法MDD


save('ca-netscience_Kendall.mat', 'ISBC_DC_K', 'ISBC_BC_K','ISBC_CC_K','ISBC_ISC_K');