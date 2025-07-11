%%
%%读取网络文件，找到最大连通分量
clear;
net=load('facebook324.txt'); %读取网络文件schematic15
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
N=size(adjMat, 1);
EWD=EWD(adjMat)';%k核分解序列值指标
burt=burt(adjMat);%形成结构洞的约束系数，越小越是结构洞
NC=NC(adjMat)';%领域度中心性
%% 计算介数中心性
BC = compute_betweenness(adjMat)';
CC=closeness(adjMat);
d = degreea( adjMat ); %度中心性
LSZ=zeros(1,N);

for i = 1:N
    LSZ(i)=log(BC(i)+5)*log(EWD(i))*exp(NC(i)/sum(d))/(1/CC(i))/burt(i);
end


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

cij=max(kcore)-min(kcore);  %KSGC算法
for i=1:N
    disp(i);
    j=1;%第一层
    nei0=[i];%第0层
    nei1=find(adjMat(i,:));%算出节点i的直接邻居
    par=[];
    nei1=find(adjMat(i,:));
    len1=length(nei1);
    par1=[];
    for j=1:len1
        par1(j)=exp((kcore(i)-kcore(nei1(j)))/cij)*d(i)*d(nei1(j));
    end

    par11=sum(par1);
    j=2;%第二层
    %     nei2=neighbor(mixedsig,nei1');
    % a=neighbor(mixedsig,nei1);
    a=[];
    for v=1:length(nei1)
        %     rr=mixedsig(nei1(v));
        k=find(adjMat(nei1(v),:));
        a=union(k,a);
    end
    c=setdiff(a,nei0);%得到节点i的二跳邻居
    nei2=setdiff(c,nei1);%得到节点i的二跳邻居

    par2=[];
    for k=1:length(nei2)
        par2(k)=exp((kcore(i)-kcore(nei2(k)))/cij)*d(i)*d(nei2(k))/2^2;
        %         par2(k)=BBB(i,nei2(k))*(1/k1)^2+BBA(i,nei2(k))*(1/k1)^3;
    end
    par22=sum(par2);

    j=3;%第3层
    a=[];
    for v=1:length(nei2)
        %     rr=mixedsig(nei1(v));
        k=find(adjMat(nei2(v),:));
        a=union(k,a);
    end
    c=setdiff(a,union(nei0,nei1));%得到节点i的3跳邻居
    nei3=setdiff(c,nei2);
    % BBA4=mixedsig^4;
    par3=[];
    for k=1:length(nei3)
        par3(k)=exp((kcore(i)-kcore(nei3(k)))/cij)*d(i)*d(nei3(k))/3^2;
    end
    par33=sum(par3);
%          runim(i)=par11+par22+par33;
   KSGC(i)=par11;
end
average_i=sir_main(adjMat,1000,0.047) ;%传播概率 每个节点的平均感染数0.047





%% 创建散点图
figure; % 创建一个新的图形窗口
scatter(d,average_i, 'r+'); 
xlabel('Degree');
ylabel('Φ', 'Rotation', 0);

figure; % 创建一个新的图形窗口
scatter(h,average_i, 'r+'); 
xlabel('H-index');
ylabel('Φ', 'Rotation', 0);

figure; % 创建一个新的图形窗口
scatter(p,average_i, 'r+'); 
xlabel('Pagerank');
ylabel('Φ', 'Rotation', 0);

figure; % 创建一个新的图形窗口
scatter(MDD_core,average_i, 'r+'); 
xlabel('MDD');
ylabel('Φ', 'Rotation', 0);

figure; % 创建一个新的图形窗口
scatter(kcore,average_i, 'r+'); 
xlabel('K-shell');
ylabel('Φ', 'Rotation', 0);

figure; % 创建一个新的图形窗口
scatter(neighbor_Core,average_i, 'r+'); 
xlabel('NCC');
ylabel('Φ', 'Rotation', 0);

figure; % 创建一个新的图形窗口
scatter(LSZ,average_i, 'r+'); 
xlabel('LSZ');
ylabel('Φ', 'Rotation', 0);

figure; % 创建一个新的图形窗口
scatter(KSGC,average_i, 'r+'); 
xlabel('KSGC');
ylabel('Φ', 'Rotation', 0);

figure; % 创建一个新的图形窗口
scatter(EWD,average_i, 'r+'); 
xlabel('NC');
ylabel('Φ', 'Rotation', 0);

