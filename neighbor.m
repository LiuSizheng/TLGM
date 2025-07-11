function a= neighbor(train, i)%train 是个邻接矩阵,函数用来返回节点i的邻居节点集合
N=size(train,1);
a=find(train(i,:));
end