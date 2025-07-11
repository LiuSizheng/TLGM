%%
%%读取网络文件，找到最大连通分量

clear;
net=load('netscience.txt'); %读取网络文件test0425
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

%%

%%计算pc值，找到初始群C，吸收节点，循环得到全部的群落Cs

A=adjMat;
N=size(A,1);  %节点数
degree=sum(A,2); %按行求和，得节点的度
avgdeg=sum(degree)/N;  %平均度
for i=1:N
    B1=find(A(i,:));%和第i个节点相连的节点索引
    len1=length(B1);%len1表示i的邻居数
    ddin(i)=0;
    for j=1:len1
        B2=find(A(B1(j),:));%i的邻居之一j的邻居
        jiaoji=intersect(B1,B2);%i与j的邻居交集
        ddin(i)=length(jiaoji)+1+ddin(i);
    end
    pc(i)=degree(i)*sum(ddin(i));
end
Cs={};
while length(find(pc==0))<N*0.6  %达到网络节点总数的60%停止

    C={};
    [a,b]=sort(pc,'descend'); %a-降序序列 b-节点索引
    [pcmax,index]=max(pc);
    biggest=find(pc==pcmax);%找到pcmax的节点
    %len2=length(biggest)
    select=degree(biggest);
    [~,max_degree]=max(select); %选择pcmax的节点中度最大的一个
    C{1}=biggest(max_degree);
    biggest_neibor=find(A(C{1},:));
    for i =1:length(biggest_neibor)
        if pc(biggest_neibor(i))>0.5*pcmax  %把pc值大于一半pcmax的加入C
            C{end+1}=biggest_neibor(i);
        end
    end
    Cpro=C;  %拓展C
    Cpromat=cell2mat(Cpro); %把元胞数组转为数组
    %avgdtmp=sum(degree(Cpromat))/length(Cpromat); %拓展C的平均度
    %num=0;
    for a=1:3
        C_neibortmp=[];
        for i =1:length(Cpromat)
            tmp=find(A(Cpromat(i),:));
            C_neibortmp=union(tmp, C_neibortmp);
        end
        C_neibor=setdiff(C_neibortmp,Cpromat); %得到C的一层邻居
        C_neibor_degree=degree(C_neibor);
        [c,d]=sort(C_neibor_degree,'descend');%按度降序
        for i = 1:length(C_neibor)
            kin=0;
            for j = 1:length(C_neibortmp)
                if A(C_neibor(d(i)),C_neibortmp(j))==1
                    kin=kin+1;
                end
            end
            kout=degree(C_neibor(d(i)))-kin;
            if kin>=kout
                Cpro{end+1}=C_neibor(d(i));
            end
            n=0;%目前Cs中已加入的节点数
            for k=1:length(Cs)
                n=n+length(Cs{k});
            end
            n=n+length(Cpro);
            if n>N*0.6
                Cpro=Cpro(1:end-1);
                break
            end
        end
        Cpromat=cell2mat(Cpro);
        if n>N*0.6
            break
        end
    end
    % while avgdtmp>=1.1*avgdeg
    %     num=num+1;
    %     C_neibortmp=[];
    %     for i =1:length(Cpromat)
    %         tmp=find(A(Cpromat(i),:));
    %         C_neibortmp=union(tmp, C_neibortmp);
    %     end
    %     C_neibor=setdiff(C_neibortmp,Cpromat); %得到C的一层邻居
    %     C_neibor_degree=degree(C_neibor);
    %     [c,d]=sort(C_neibor_degree,'descend');%按度降序
    %     for i = 1:length(C_neibor)
    %         kin=0;
    %         for j = 1:length(C_neibortmp)
    %             if A(C_neibor(d(i)),C_neibortmp(j))==1
    %                 kin=kin+1;
    %             end
    %         end
    %         kout=degree(C_neibor(d(i)))-kin;
    %         if kin>=kout
    %             Cpro{end+1}=C_neibor(d(i));
    %         end
    %         Cpromat=cell2mat(Cpro); %把元胞数组转为数组
    %         avgdtmp=sum(degree(Cpromat))/length(Cpromat);
    %         if avgdtmp<1.1*avgdeg   %当添加节点后平均度小于规定，就停止，且移除
    %             Cpro=Cpro(1:end-1);
    %             break
    %         end
    %         n=0;%目前Cs中已加入的节点数
    %         for k=1:length(Cs)
    %             n=n+length(Cs{k});
    %         end
    %         n=n+length(Cpro);
    %         if n>N*0.6
    %             Cpro=Cpro(1:end-1);
    %             break
    %         end
    %     end
    %     if n>N*0.6
    %         break
    %     end
    % end
    Cpromat=cell2mat(Cpro);
    pc(Cpromat)=0;   %已经聚类的节点pc置为零
    non_zero_rows = find(pc ~= 0);
    zero_cols = find(pc == 0);
    for i = 1:length(non_zero_rows)  %邻接矩阵断掉已经聚类节点的连接，方便循环中邻居的判定
        for j = 1:length(zero_cols)
            A(non_zero_rows(i), zero_cols(j)) = 0;
            A(zero_cols(j), non_zero_rows(i)) = 0;
        end
    end
    Cs{end+1}=Cpro;
end

%%

%%找到群落中度最大的节点

C=adjMat;
degree1=sum(C,2); %按行求和，得节点的度
for i=1:length(Cs)
    [dmax,didx]=max(degree1(cell2mat(Cs{i})));
    dmaxnode(i)=Cs{i}{didx};%找到每个簇中度最大的节点
end

%%

%%节点局部影响计算
kcore=core_numbers(sparse(adjMat))';
LI=[];
for i=1:N
    B3=find(C(i,:));%和第i个节点相连的节点索引
    len2=length(B3);%len3表示i的邻居数
    ks=0;
    for j=1:len2
        ks=ks+kcore(B3(j));
    end
    LI(end+1)=ks;
end
shortest_distances = dijkstra_for_targets(adjMat,dmaxnode);%找到全局关键节点到其他节点的最短距离

%%

%%总体影响力的计算
for i=1:N
    a=0;
    for j =1:length(dmaxnode)
        a=a+LI(dmaxnode(j))/(2^(shortest_distances(j,i)));
    end
    GAC(i)=LI(i)*a;
end




%%


d = degreea( adjMat ); %度中心性

h = h_index(adjMat); %H指数

p=pagerank(adjMat,N);%PageRank

MDD_core=MDD(adjMat);%混合度分解法（MDD）

kcore=core_numbers(sparse(adjMat))'; %k-壳分解

[r,l]=size(adjMat);  %邻域核度NCC
neighbor_Core=[];
for i=1:r
    B=find(adjMat(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_Core(i)=0;
    for j=1:len1
        neighbor_Core(i)=kcore(B(j))+neighbor_Core(i);
    end
end

average_i=sir_main(adjMat,1000,0.125) ;%传播概率 每个节点的平均感染数

save('netscience.mat', 'd', 'h','p','MDD_core','kcore','neighbor_Core','GAC','average_i');




