function c=pagerank(A,n)   %传入邻接图A和网络规模n
%1.计算每个点的出度
outDegree=zeros(n,1);
for i=1:n
    current_outDegree=0;
    for j=1:n
        if A(i,j)>0
            current_outDegree=current_outDegree+1;
        end
    end
    outDegree(i,1)=current_outDegree;
end
 
%2.基本的pagerank算法，得到的google矩阵A1如下
A1=A;
for i=1:n
    for j=1:n
        if outDegree(i,1)==0                 %如果一个点出度为0，为了防止这个点把整个网路的PR值耗散掉，设定这个点没有链接的情况下，随机跳转到其他网络
            A1(i,j)=1/n;
        else
        A1(i,j)=A1(i,j)/outDegree(i,1);
        end
    end
end
A1
%3.初始化PR向量，使得浏览每个页面的概率相同
PR=zeros(n,1);
for i=1:n
    PR(i,1)=1/n;
end
PR
%4.把谷歌矩阵A1转置，右边乘上PR矩阵，就是一次迭代
%迭代收敛的意思是，这次迭代和下次迭代结果相同，就停止
num=0;
runpagerank(A1,PR)
while 1
    if round(PR,4)==round(runpagerank(A1,PR),4) %这里如果不用估值（前4），一直没法终止算法
        c=PR;
        break;
    else
        PR=runpagerank(A1,PR);
        PR
        num=num+1
    end
end
