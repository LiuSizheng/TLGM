%%
%%读取网络文件，找到最大连通分量

clear;
net=load('protain.txt'); %读取网络文件test0425
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

mix=adjMat;
N=size(mix,1);  %节点数
degree=sum(mix,2)'; %按行求和，得节点的度
avgdeg=sum(degree)/N;  %平均度
for i=1:N
    B1=find(mix(i,:));%和第i个节点相连的节点索引
    len1=length(B1);%len1表示i的邻居数
    ddin(i)=0;
    for j=1:len1
        B2=find(mix(B1(j),:));%i的邻居之一j的邻居
        jiaoji=intersect(B1,B2);%i与j的邻居交集
        ddin(i)=length(jiaoji)+ddin(i);
    end
    pc(i)=degree(i)*sum(ddin(i));
end
core=pc;
C={};
ii=1;
totlenc=0;
num=0;
threshold=0;
while threshold<=0.6  %达到网络节点总数的60%停止

    
    %[a,b]=sort(core,'descend'); %a-降序序列 b-节点索引
    [pcmax,index]=max(core);
    biggest=find(core==pcmax);%找到pcmax的节点
    aaa=biggest(find(degree(biggest)==max(degree(biggest))));%找到pc值最大的里面度最大的
    len0=length(aaa);
    initial_set(ii)=aaa(ceil(rand(1)*len0));%随机取得最大pc值中的一个作为起始节点%%%%
    num=num+1;
    if num/N<=0.6
        %initial_set(ii)=aaa(1);
        C{ii}=initial_set(ii);
        neighbor_of_i=find(mix(initial_set(ii),:))%初始节点的邻居都找出来
        zz=find(pc(neighbor_of_i)>=pc(initial_set(ii))/2);%找到pc值大于1/2pcmax的节点
        pc_of_zz=pc(neighbor_of_i(zz));
        togther0=[pc_of_zz;neighbor_of_i(zz)]';
        togtherr0=sortrows(togther0)';
        for j=length(zz):-1:1
            C{ii}=[C{ii},togtherr0(2,j)];
            num=num+1;
            if num/N>0.6
                C{ii} = C{ii}(1:end-1);
                break
            end
        end

        % C{ii}=neighbor_of_i(zz);%找到节点i的内部链接的各个节点,组成子图C
        % num=num+length(C{ii});
        %C{ii}=[C{ii},initial_set(ii)];
        if num/N<=0.6
            first_neighbor_of_C=neighbor_of_node_set(mix,C{ii});
            degree_of_first_neighbor_of_C=degree(first_neighbor_of_C);
            togther=[degree_of_first_neighbor_of_C;first_neighbor_of_C]';
            togtherr=sortrows(togther)';
            
            %for j=length(first_neighbor_of_C):-1:1
            for j=1:length(first_neighbor_of_C)
                ddd=neighbor(mix,togtherr(2,j));
                jiaoji=intersect(ddd,C{ii});
                chaji=setdiff(ddd,jiaoji);
                if length(jiaoji)>length(chaji)
                    C{ii}=[C{ii},togtherr(2,j)];
                    num=num+1;
                    if num/N>0.6
                        C{ii} = C{ii}(1:end-1);
                        break
                    end
                end
            end
        else
            break
        end
        if num/N<=0.6
            second_neighbor_of_C=neighbor_of_node_set(mix,C{ii});
            degree_of_second_neighbor_of_C=degree(second_neighbor_of_C);
            togther=[degree_of_second_neighbor_of_C;second_neighbor_of_C]';
            togtherr=sortrows(togther)';
            %for j=length(second_neighbor_of_C):-1:1
            for j=1:length(second_neighbor_of_C)
                ddd=neighbor(mix,togtherr(2,j));
                jiaoji=intersect(ddd,C{ii});
                chaji=setdiff(ddd,jiaoji);
                if length(jiaoji)>length(chaji)
                    C{ii}=[C{ii},togtherr(2,j)];
                    num=num+1;
                    if num/N>0.6
                        C{ii} = C{ii}(1:end-1);
                        break
                    end
                end
        
            end
        else
            break
        end
        if num/N<=0.6
            third_neighbor_of_C=neighbor_of_node_set(mix,C{ii});
            degree_of_third_neighbor_of_C=degree(third_neighbor_of_C);
            togther=[degree_of_third_neighbor_of_C;third_neighbor_of_C]';
            togtherr=sortrows(togther)';
            %for j=length(third_neighbor_of_C):-1:1
            for j=1:length(third_neighbor_of_C)
                ddd=neighbor(mix,togtherr(2,j));
                jiaoji=intersect(ddd,C{ii});
                chaji=setdiff(ddd,jiaoji);
                if length(jiaoji)>length(chaji)
                    C{ii}=[C{ii},togtherr(2,j)];
                    num=num+1;
                    if num/N>0.6
                        C{ii} = C{ii}(1:end-1);
                        break
                    end
                end
        
            end
        else
            break
        end
    
        %C_number_of_initial_i(i)=length(C{i});
        linshi3=C{ii};
        lenc=length(linshi3);
        totlenc=totlenc+lenc;
        linshi4=linshi3(find(degree(linshi3)==max(degree(linshi3))));
        initial_set(ii)=linshi4(ceil(rand(1)*length(linshi4)));
        %initial_set(ii)=linshi4(1);
    
        mix(C{ii},:)=0;
        mix(:,C{ii})=0;
        degree=sum(mix);
        MN=sparse(mix);
        for i=1:N
            B1=find(mix(i,:));%和第i个节点相连的节点索引
            len1=length(B1);%len1表示i的邻居数
            ddin(i)=0;
            for j=1:len1
                B2=find(mix(B1(j),:));%i的邻居之一j的邻居
                jiaoji=intersect(B1,B2);%i与j的邻居交集
                ddin(i)=length(jiaoji)+ddin(i);
            end
            core(i)=degree(i)*sum(ddin(i));
        end
        
        
        ii=ii+1;
        threshold=totlenc/N;    
    else
        initial_set=initial_set(1:end-1);
        break
    end
end



%%节点局部影响计算
kcore=core_numbers(sparse(adjMat))';
LI=[];
for i=1:N
    B3=find(adjMat(i,:));%和第i个节点相连的节点索引
    len2=length(B3);%len3表示i的邻居数
    ks=0;
    for j=1:len2
        ks=ks+kcore(B3(j));
    end
    LI(end+1)=ks;
end
shortest_distances = dijkstra_for_targets(adjMat,initial_set);%找到全局关键节点到其他节点的最短距离

%%

sd=shortest_distances';
for i=1:N
    GAC(i)=0;
    for j=1:length(initial_set)
        GAC(i)=LI(i)*LI(initial_set(j))/(2^sd(i,j))+GAC(i);
    end
end




%%

%totalNumbers = sum(cellfun(@length, C)); % 


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
average_i=sir_main(adjMat,1000,0.063) ;%传播概率 每个节点的平均感染数0.047

save('protain.mat', 'd', 'h','p','MDD_core','kcore','neighbor_Core','GAC','KSGC','average_i');

average_iall={};
w = waitbar(0, 'Please wait...'); % 初始化进度条

for i = 1:15
    waitbar(i/15,w , sprintf('Progress: %d %%', floor(i/15*100))); % 更新进度条
    disp(i);
    average_iall{i} = sir_main(adjMat, 1000, 0.01*i);
end

close(w); % 运行结束后关闭进度条

%a00_hefei1=[];a00_h=[];a00_pagerank=[];a00_M_core=[];a00_neighbor_core=[];a00_degree=[];a0_core=[];a00_KSGC=[];a00_h=[];
for i=1:15
    a00_degree(i)=corr(average_iall{i}',d','type','Kendall');%度中心性
    a00_h_index(i)=corr(average_iall{i}',h,'type','Kendall');%H指数
    a00_pagerank(i)=corr(average_iall{i}',p,'type','Kendall');%PageRank中心性
    a00_MDD_core(i)=corr(average_iall{i}',MDD_core','type','Kendall');%混合度分解法MDD
    a00_GAC(i)=corr(average_iall{i}',GAC','type','Kendall');
    a00_kcore(i)=corr(average_iall{i}',kcore','type','Kendall');%k壳分解法
    a00_neighbor_core(i)=corr(average_iall{i}',neighbor_Core','type','Kendall');%邻域核心度中心性（NCC）
    a00_KSGC(i)=corr(average_iall{i}',KSGC','type','Kendall');%基于 K 核分解法的引力模型改进方法KSGC
end

save('protain_sir1_15_Kendall.mat', 'a00_degree', 'a00_h_index','a00_pagerank','a00_MDD_core','a00_GAC','a00_kcore','a00_neighbor_core','a00_KSGC');






