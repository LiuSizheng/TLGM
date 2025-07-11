clear ;
A=load('netscience.txt');
%A=A+1;
r=0.2;
%% for linjiejuzhen
TT=A(:, 1:2);
mixedsig=zeros(max(max(TT)));
N=size(mixedsig,1);
len=length(TT);
for i=1:len
    mixedsig(TT(i,1),TT(i,2))=1;
    mixedsig(TT(i,2),TT(i,1))=1;
end
%%%for linjiejuzhen

%%for largestcomponent
kk=sparse(mixedsig);
[a b]=components(kk);
[B]=largestcomponent(mixedsig);
mixedsig=mixedsig(B,B);
%%for largestcomponent


 A=mixedsig;
  N=size(A,1); 
  b=sum(A,2);
  k1=sum(b)/N;
  k2=sum(b'*b)/N; 
irate=k1/k2;
core=core_numbers(sparse(mixedsig))';
ccfs = sum(clustering_coefficients(sparse(mixedsig)))/N;
dmean=sum(degreea( mixedsig))/N;
[D,aver_D]=Aver_Path_Length(mixedsig);
d = degreea( mixedsig );
k=sum(d)/N;
mix=mixedsig;
M_core=MDD(mixedsig);
%bet=betweenness_centrality(sparse(mixedsig));
%clo=Closeness_Centrality(mixedsig);
%result = find_shortest_distances(N,len, mixedsig);
for i=1:N
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    ddin(i)=0;
    for j=1:len1
         B2=find(mixedsig(B(j),:));
          jiaoji=intersect(B,B2);
        ddin(i)=d(i)*length(jiaoji)+ddin(i);
    end
end


C={};
[a,b]=sort(ddin,'descend');
[max_core,indexm]=max(ddin);%找出第i轮中得分最高的那个值
biggest_core=find(ddin==max_core);
len0=length(biggest_core);
DegreeAll=sum(mixedsig);
m=mixedsig;
[r,l]=size(mixedsig);
neighbor_Core=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_Core(i)=0;
    for j=1:len1
        neighbor_Core(i)=core(B(j))+neighbor_Core(i);
    end
end
 [r,l]=size(mixedsig);
extend_neighbor_core=[];
 for i=1:r
 B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
   extend_neighbor_core(i)=0;
     for j=1:len1
    extend_neighbor_core(i)=neighbor_Core(B(j))+extend_neighbor_core(i);
     end
 end
 
 CC=burt(mixedsig);
CC=-1*CC;
CCS=burt_spreading_out(mixedsig);
 for i=1:N
     %extend_cluster_coreness(i)=core(i)+neighbor_core(i)/power(1.5,C(i));
     coreness_CCS(i)=core(i)/power(10,CCS(i));
 end
k=30;
K=30;  %%%节点的初始激活数目为K
threshold=0;
totlenc=0;
i=1;
while threshold<=0.6
%for i=1:K
    aaa=biggest_core(find(DegreeAll(biggest_core)==max(DegreeAll(biggest_core))));
    len0=length(aaa);
    initial_set(i)=aaa(ceil(rand(1)*len0));%随机取得最大核数中的一个作为起始节点%%%%
    neighbor_of_i=find(mix(initial_set(i),:))%初始节点的邻居都找出来
    zz=find(ddin(neighbor_of_i)>=ddin(initial_set(i))/2);
    C{i}=neighbor_of_i(zz);%找到节点i的内部链接的各个节点,组成子图C
    C{i}=[C{i},initial_set(i)];
    first_neighbor_of_C=neighbor_of_node_set(mix,C{i});
    degree_of_first_neighbor_of_C=DegreeAll(first_neighbor_of_C);
    togther=[degree_of_first_neighbor_of_C;first_neighbor_of_C]';
    togtherr=sortrows(togther)';
    for j=1:length(first_neighbor_of_C)
        ddd=neighbor(mix,togtherr(2,j));
        jiaoji=intersect(ddd,C{i});
        chaji=setdiff(ddd,jiaoji);
        if length(jiaoji)>length(chaji)
            C{i}=[C{i},togtherr(2,j)];
        end
    end
    
    second_neighbor_of_C=neighbor_of_node_set(mix,C{i});
    degree_of_second_neighbor_of_C=DegreeAll(second_neighbor_of_C);
    togther=[degree_of_second_neighbor_of_C;second_neighbor_of_C]';
    togtherr=sortrows(togther)';
    for j=1:length(second_neighbor_of_C)
        ddd=neighbor(mix,togtherr(2,j));
        jiaoji=intersect(ddd,C{i});
        chaji=setdiff(ddd,jiaoji);
        if length(jiaoji)>length(chaji)
            C{i}=[C{i},togtherr(2,j)];
        end
        
    end
    third_neighbor_of_C=neighbor_of_node_set(mix,C{i});
    degree_of_third_neighbor_of_C=DegreeAll(third_neighbor_of_C);
    togther=[degree_of_third_neighbor_of_C;third_neighbor_of_C]';
    togtherr=sortrows(togther)';
    for j=1:length(third_neighbor_of_C)
        ddd=neighbor(mix,togtherr(2,j));
        jiaoji=intersect(ddd,C{i});
        chaji=setdiff(ddd,jiaoji);
        if length(jiaoji)>length(chaji)
            C{i}=[C{i},togtherr(2,j)];
        end
        
    end
    C_number_of_initial_i(i)=length(C{i});
    linshi3=C{i};
    lenc=length(linshi3);
    totlenc=totlenc+lenc;
    linshi4=linshi3(find(DegreeAll(linshi3)==max(DegreeAll(linshi3))));
    initial_set(i)=linshi4(ceil(rand(1)*length(linshi4)));
    
    mix(C{i},:)=0;
    mix(:,C{i})=0;
    [a,b]=sort(ddin,'descend');
    MN=sparse(mix);
    core=iniddin(mix,N);
    %core=core_numbers(MN)';
    [a,b]=sort(core,'descend');
    [max_core,indexm]=max(core);%找出第i轮中得分最高的那个值
    biggest_core=find(core==max_core);
    len0=length(biggest_core);
    DegreeAll=sum(mix);
    i=i+1;  
    threshold=totlenc/N;
end
lenghefei=length(initial_set);
core=core_numbers(sparse(mixedsig))';
 h = (h_index(mixedsig))';
 neighbor_h=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_h(i)=0;
    for j=1:len1
        neighbor_h(i)=core(B(j))+neighbor_h(i);
    end
end
shortest_distances = dijkstra_for_targets(mixedsig, initial_set);
sd=shortest_distances';
for i=1:N
    hefei1(i)=0;
    for j=1:lenghefei
    %DD = dijkstra(mixedsig,initial_set(j));
    %sssss=core(i);
    %dfs3e=core(initial_set(j));
    %    hefei1(i)=neighbor_Core(i)*neighbor_Core(initial_set(j))/(2^result(i,initial_set(j)))+hefei1(i);
            hefei1(i)=neighbor_Core(i)*neighbor_Core(initial_set(j))/(2^sd(i,j))+hefei1(i);
                    %    hefei1(i)=(neighbor_h(i)+neighbor_Core(i))/1*(neighbor_h(initial_set(j))*neighbor_Core(initial_set(j))/1/(2^sd(i,j))+hefei1(i));

               % hefei1(i)=exp(neighbor_Core(i))*exp(neighbor_Core(initial_set(j)))/(2^sd(i,j))+hefei1(i);

    end
end

 [r,l]=size(mixedsig);
extend_hefei=[];
 for i=1:r
 B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
   extend_hefei(i)=0;
     for j=1:len1
    extend_hefei(i)=hefei1(B(j))+extend_hefei(i);
     end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    disp(i);
    j=1;%第一层
    nei0=[i];%第0层
    nei1=find(mixedsig(i,:));%算出节点i的直接邻居
    par=[];
    nei1=find(mixedsig(i,:));
    len1=length(nei1);
    par1=[];
    for j=1:len1
        par1(j)=core(i)*core(nei1(j));
    end
    
    par11=sum(par1);
    j=2;%第二层
    %     nei2=neighbor(mixedsig,nei1');
    % a=neighbor(mixedsig,nei1);
    a=[];
    for v=1:length(nei1)
        %     rr=mixedsig(nei1(v));
        k=find(mixedsig(nei1(v),:));
        a=union(k,a);
    end
    c=setdiff(a,nei0);%得到节点i的二跳邻居
    nei2=setdiff(c,nei1);%得到节点i的二跳邻居
    
    par2=[];
    for k=1:length(nei2)
        par2(k)=core(i)*core(nei2(k))/2^2;
        %         par2(k)=BBB(i,nei2(k))*(1/k1)^2+BBA(i,nei2(k))*(1/k1)^3;
    end
    par22=sum(par2);
    
    j=3;%第3层
    a=[];
    for v=1:length(nei2)
        %     rr=mixedsig(nei1(v));
        k=find(mixedsig(nei2(v),:));
        a=union(k,a);
    end
    c=setdiff(a,union(nei0,nei1));%得到节点i的3跳邻居
    nei3=setdiff(c,nei2);
    % BBA4=mixedsig^4;
    par3=[];
    for k=1:length(nei3)
        par3(k)=core(i)*core(nei3(k))/3^2;
    end
    par33=sum(par3);
    runk(i)=par11+par22+par33;
end
cij=max(core)-min(core);
for i=1:N
    disp(i);
    j=1;%第一层
    nei0=[i];%第0层
    nei1=find(mixedsig(i,:));%算出节点i的直接邻居
    par=[];
    nei1=find(mixedsig(i,:));
    len1=length(nei1);
    par1=[];
    for j=1:len1
        par1(j)=exp((core(i)-core(nei1(j)))/cij)*d(i)*d(nei1(j));
    end
    
    par11=sum(par1);
    j=2;%第二层
    %     nei2=neighbor(mixedsig,nei1');
    % a=neighbor(mixedsig,nei1);
    a=[];
    for v=1:length(nei1)
        %     rr=mixedsig(nei1(v));
        k=find(mixedsig(nei1(v),:));
        a=union(k,a);
    end
    c=setdiff(a,nei0);%得到节点i的二跳邻居
    nei2=setdiff(c,nei1);%得到节点i的二跳邻居
    
    par2=[];
    for k=1:length(nei2)
        par2(k)=exp((core(i)-core(nei2(k)))/cij)*d(i)*d(nei2(k))/2^2;
        %         par2(k)=BBB(i,nei2(k))*(1/k1)^2+BBA(i,nei2(k))*(1/k1)^3;
    end
    par22=sum(par2);
    
    j=3;%第3层
    a=[];
    for v=1:length(nei2)
        %     rr=mixedsig(nei1(v));
        k=find(mixedsig(nei2(v),:));
        a=union(k,a);
    end
    c=setdiff(a,union(nei0,nei1));%得到节点i的3跳邻居
    nei3=setdiff(c,nei2);
    % BBA4=mixedsig^4;
    par3=[];
    for k=1:length(nei3)
        par3(k)=exp((core(i)-core(nei3(k)))/cij)*d(i)*d(nei3(k))/3^2;
    end
    par33=sum(par3);
%          runim(i)=par11+par22+par33;
   KSGC(i)=par11;
end
h = h_index(mixedsig);
c=pagerank(mixedsig,N);
average_i=sir_main(mixedsig,1000,0.063) ;
%%%%%%%%%%%%%%%%%%%%%
average_iall={};
for i=1:26
    disp(i);
    average_iall{i}=sir_main(mixedsig,1000,0.01*i);
end
a00_hefei1=[];a00_h=[];a00_pagerank=[];a00_M_core=[];a00_neighbor_core=[];a00_degree=[];a0_core=[];a00_KSGC=[];a00_h=[];
for i=1:15
%average_i=sir_main(mixedsig,1000,0.05) ;
%average_i=average_iall(1) ;
a00_degree(i)=corr(average_iall{i}',d','type','Kendall');
a00_h(i)=corr(average_iall{i}',h,'type','Kendall');
a00_pagerank(i)=corr(average_iall{i}',c,'type','Kendall');
a00_M_core(i)=corr(average_iall{i}',M_core','type','Kendall');
a00_bet(i)=corr(average_iall{i}',bet,'type','Kendall');
a00_clo(i)=corr(average_iall{i}',clo','type','Kendall');
%a00_densityim=corr(average_i',runkim','type','Kendall');
a00_hefei1(i)=corr(average_iall{i}',hefei1','type','Kendall');
a00_extend_hefei(i)=corr(average_iall{i}',extend_hefei','type','Kendall');
a0_core(i)=corr(average_iall{i}',core','type','Kendall');
a00_extend_neighbor_core(i)=corr(average_iall{i}',extend_neighbor_core','type','Kendall')
a00_neighbor_core(i)=corr(average_iall{i}',neighbor_Core','type','Kendall');
a00_runk(i)=corr(average_iall{i}',runk','type','Kendall');
a00_KSGC(i)=corr(average_iall{i}',KSGC','type','Kendall');
end
all_average_sir=[a00_hefei1;a00_M_core;a00_bet;a00_clo;a00_extend_hefei;a00_neighbor_core;a00_degree;a0_core;a00_extend_neighbor_core;a00_runk;a00_KSGC];
all_average_sir2=[a00_hefei1;a00_h;a00_pagerank;a00_M_core;a00_neighbor_core;a00_degree;a0_core;a00_KSGC];
all_average_sir3=all_average_sir2';
save('LFR2000-k10_average_sir1_15.mat');
all_average_sir3=all_average_sir2';

