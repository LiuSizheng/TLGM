 
clear all;
%irate=0.15;
%A=load('erdos.txt');
A=load('URVemail.txt');
%A=load('blog.txt');


 TT=A(:, 1:2);
 mixedsig=zeros(max(max(TT)));
 len=length(TT);
 for i=1:len
     mixedsig(TT(i,1),TT(i,2))=1;
     mixedsig(TT(i,2),TT(i,1))=1;
 end
kk=sparse(mixedsig);
[a b]=components(kk);
 [B] = largestcomponent(mixedsig);
%kk=sparse(mixedsig);
%[a b]=components(kk);
% [B] = largestcomponent(mixedsig);
mixedsig=mixedsig(B,B);
 A=mixedsig;
  N=size(A,1); 
  b=sum(A,2);
  k1=sum(b)/N;
  k2=sum(b'*b)/N; 
irate=k1/k2;
average_i=sir_main(mixedsig,500) ;
 %先求出各个节点的集聚系数
 N=size(mixedsig,2);
C=zeros(1,N);
for i=1:N
    a=find(mixedsig(i,:)==1); %寻找子图的邻居节点
    b=find(mixedsig(:,i)==1);
    m=union(a,b'); 
    k=length(m);
    if k==1
        disp(['节点',int2str(i),'只有一个邻居节点，其聚类系数为0']);
        C(i)=0;
    else
        B=mixedsig(m,m);
        C(i)=length(find(B==1))/(k*(k-1));
    end
end
aver_C=mean(C);
max_C=max(C);


%average_i=[5.29000000000000	4.88000000000000	4.780000000000000	4.53000000000000	3.90000000000000	3.79000000000000	3.05000000000000	3.04000000000000	2.50000000000000	2.51000000000000	2.09000000000000	2.05000000000000	2.02000000000000	2.34000000000000	1.61000000000000	4.34000000000000	2.19000000000000	2.22000000000000	2.21000000000000	2.24000000000000	2.23000000000000	2.20000000000000	2.67000000000000	2.32000000000000	1.56000000000000	1.55000000000000
%];
j=1;   %第j层
t=1;  %
r=1;
tad=mixedsig;  %邻接矩阵
k=0;
ii=1;
ret_m=[];  %返回值；  i，j  ： 第i层：所有节点
ret_r=[]; 
len=length(mixedsig);  %矩阵长度
tt=1;
n=1;
nn=1;
m_t=[];
m_tt=[];
m_ttt=[];
p=1;
while (tt==1)   %控制最高层结束
    disp(tt);
    if sum(sum(tad))==0  %所有元素为0，则退出
        break;
    end
    t=1;  %控制第j层  计算
    while(t==1)  %  每循环一次，tad改变，去掉度数小于j的节点；直到没有度数小于j的节点
        t=0;
        ii=1; %第j层第ii个节点
        for i=1:len  %从矩阵1至len行,  去掉小于j的节点
            k=sum(tad(i,:));   %计算i行度数
            if k==0  %度数为0，下一个i值
                % t=1 ;  %
                continue;
            elseif k<=j  %度数小于j层
                t=1; %控制下一次还要循环
                tad(i,:)=0;  %i节点加到j层，将度数至为0，所有i行值至为0
                tad(:,i)=0; %相应i列至为0

                ret_m(j,ii)=i;   %将i节点加到j层
               ret_r(j,r)=i;
               r=r+1;
                m_ttt=union(ret_m(j,ii),m_ttt);
                for n=1:len   %判断其他行有无因为 将tad(:,i)元素设为0 而所有行变0
                    if sum(tad(n,:))==0
                        m_t=n;
                        m_tt = intersect(m_ttt,m_t);
                        if  length(m_tt)==0  %length(m_t)~=0 &&
                            ii=ii+1;
                            r=r+1;
                            ret_m(j,ii)=n;
                            ret_r(j,r)=n;
                            m_ttt=union(ret_m(j,ii),m_ttt);
                            m_t=[];
                             disp(r);
                        end
                    end



                end %for n=1:len
            end % if k==0
             ii=ii+1;
             r=r+1;
        end %end of  i=1:len
         
    end  % end of while(t==1)
    disp(ret_r);
    j=j+1;
    ii=1;
    r=1;
end   % end of while(tt==1)
%disp(ret_m);


%接下来将核数保留在数组core里
core=[];
[r,l]=size(ret_r);
for i=1:r
B=find(ret_r(i,:));
k1=length(B);
for j=1:k1
    core(ret_r(i,B(j)))=i;
end
end

 
 
 %改进算法coreness_burt
 N=size(mixedsig,2);
 cc=burt(mixedsig);
 for i=1:N
     %extend_cluster_coreness(i)=core(i)+neighbor_core(i)/power(1.5,C(i));
     coreness_burt(i)=core(i)/power(10,cc(i));
 end
 %编写coreness-burt的扩展算法
 [r,l]=size(mixedsig);
extend_coreness_burt=[];
 for i=1:r
 B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
   extend_coreness_burt(i)=0;
     for j=1:len1  
    extend_coreness_burt(i)=coreness_burt(B(j))+extend_coreness_burt(i);
     end
 end
 %下一步计算出每一个节点的邻居核数neighbor_core
m=mixedsig;
[r,l]=size(mixedsig);
neighbor_core=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_core(i)=0;
     for j=1:len1  
    neighbor_core(i)=core(B(j))+neighbor_core(i);
    end
    
end
%接下来编写扩展算法extend_neighbor_core
 [r,l]=size(mixedsig);
extend_neighbor_core=[];
 for i=1:r
 B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
   extend_neighbor_core(i)=0;
     for j=1:len1  
    extend_neighbor_core(i)=neighbor_core(B(j))+extend_neighbor_core(i);
     end
 end 


%下一步计算出改进后的每一个节点的邻居核数neighbor_burt_core
m=mixedsig;
[r,l]=size(mixedsig);
neighbor_burt_core=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_burt_core(i)=0;
     for j=1:len1  
    neighbor_burt_core(i)=coreness_burt(B(j))+neighbor_burt_core(i);
    end
    
end
%改进后的扩展算法extend_neighbor_burt_core
 [r,l]=size(mixedsig);
extend_neighbor_burt_core=[];
 for i=1:r
 B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
   extend_neighbor_burt_core(i)=0;
     for j=1:len1  
    extend_neighbor_burt_core(i)=neighbor_burt_core(B(j))+extend_neighbor_burt_core(i);
     end
 end
qqqq0_core=corr(average_i',core','type','Kendall');
qqqq1_core_burt=corr(average_i',coreness_burt','type','Kendall');
qqqq2_neighbor_burt_core=corr(average_i',neighbor_burt_core','type','Kendall');
qqqq3_extend_neighbor_burt_core=corr(average_i',extend_neighbor_burt_core','type','Kendall');
qqqq4_neighbor_core=corr(average_i',neighbor_core','type','Kendall');
qqqq5_extend_neighbor_core=corr(average_i',extend_neighbor_core','type','Kendall');
 