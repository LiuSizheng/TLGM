
clear ;
A=load('netscience.txt');
A=A+1;
r=0.2;
TT=A(:, 1:2);
mixedsig=zeros(max(max(TT)));
N=size(mixedsig,1);
len=length(TT);
for i=1:len
    mixedsig(TT(i,1),TT(i,2))=1;
    mixedsig(TT(i,2),TT(i,1))=1;
end
kk=sparse(mixedsig);
[a b]=components(kk);
[B]=largestcomponent(mixedsig);
mixedsig=mixedsig(B,B);

 A=mixedsig;
  N=size(A,1); 
  b=sum(A,2);
  k1=sum(b)/N;
  k2=sum(b'*b)/N; 
irate=k1/k2;
core=core_numbers(sparse(mixedsig))';
d = degreea( mixedsig );
k=sum(d)/N;
cdt=core_discount_two(mixedsig,core,k)
[r,l]=size(mixedsig);
neighbor_cdt=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_cdt(i)=0;
    for j=1:len1
        neighbor_cdt(i)=cdt(B(j))+neighbor_cdt(i);
    end
end


%cc=newburt(mixedsig,core);
CC=burt(mixedsig);
CC=-1*CC;
CCS=burt_spreading_out(mixedsig);
%下一步计算出每一个节点的邻居核数neighbor_CCS
m=mixedsig;
[r,l]=size(mixedsig);
neighbor_CCS=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_CCS(i)=0;
    for j=1:len1
        neighbor_CCS(i)=CCS(B(j))+neighbor_CCS(i);
    end
end

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


A=mixedsig;
N=size(A,1);
b=sum(A,2);
k1=sum(b)/N;
k2=sum(b'*b)/N;
irate=k1/k2;



%mdd=MDD(mixedsig);
%mddi=MDDI(mixedsig);
%kl=kliu(mixedsig,core);


 %改进算法coreness_burt
 N=size(mixedsig,2);

 for i=1:N
     %extend_cluster_coreness(i)=core(i)+neighbor_core(i)/power(1.5,C(i));
     coreness_CCS(i)=core(i)/power(10,CCS(i));
 end
 
%改进算法coreness_newburt
 N=size(mixedsig,2);

 for i=1:N
     %extend_cluster_coreness(i)=core(i)+neighbor_core(i)/power(1.5,C(i));
     coreness_newburt(i)=core(i)/power(10,CC(i));
 end
%  %编写coreness-burt的扩展算法
%  [r,l]=size(mixedsig);
% extend_coreness_burt=[];
%  for i=1:r
%  B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
%     len1=length(B);%len1表示数组B的长度
%    extend_coreness_burt(i)=0;
%      for j=1:len1
%     extend_coreness_burt(i)=coreness_burt(B(j))+extend_coreness_burt(i);
%      end
%  end
%  %下一步计算出每一个节点的邻居核数neighbor_core
% m=mixedsig;
% [r,l]=size(mixedsig);
% neighbor_core=[];
% for i=1:r
%     B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
%     len1=length(B);%len1表示数组B的长度
%     neighbor_core(i)=0;
%      for j=1:len1
%     neighbor_core(i)=core(B(j))+neighbor_core(i);
%     end
%
% end
%
%
%接下来编写扩展算法extend_neighbor_core
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

 
 %接下来编写扩展算法extend_neighbor_core
 [r,l]=size(mixedsig);
extend_neighbor_cdt=[];
 for i=1:r
 B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
   extend_neighbor_cdt(i)=0;
     for j=1:len1
    extend_neighbor_cdt(i)=neighbor_cdt(B(j))+extend_neighbor_cdt(i);
     end
 end
%
% %下一步计算出改进后的每一个节点的邻居核数neighbor_burt_core
% m=mixedsig;
% [r,l]=size(mixedsig);
% neighbor_burt_core=[];
% for i=1:r
%     B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
%     len1=length(B);%len1表示数组B的长度
%     neighbor_burt_core(i)=0;
%      for j=1:len1
%     neighbor_burt_core(i)=coreness_burt(B(j))+neighbor_burt_core(i);
%     end
%
% end
% %改进后的扩展算法extend_neighbor_burt_core
%  [r,l]=size(mixedsig);
% extend_neighbor_burt_core=[];
%  for i=1:r
%  B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
%     len1=length(B);%len1表示数组B的长度
%    extend_neighbor_burt_core(i)=0;
%      for j=1:len1
%     extend_neighbor_burt_core(i)=neighbor_burt_core(B(j))+extend_neighbor_burt_core(i);
%      end
%  end

average_i=sir_main(mixedsig,500,0.1) ;


a00_cdt=corr(average_i',cdt','type','Kendall');
a00_neighbor_cdt=corr(average_i',neighbor_cdt','type','Kendall');
a00_extend_neighbor_core=corr(average_i',extend_neighbor_core','type','Kendall')
a00_extend_neighbor_cdt=corr(average_i',extend_neighbor_cdt','type','Kendall')
a00_neighbor_core=corr(average_i',neighbor_Core','type','Kendall');
a00_degree=corr(average_i',d','type','Kendall');
a0_core=corr(average_i',core','type','Kendall');
a0_CC=corr(average_i',CC','type','Kendall');
a0_CCS=corr(average_i',CCS','type','Kendall');
a0_neighbor_CCS=corr(average_i',neighbor_CCS','type','Kendall');
a0_coreness_CCS=corr(average_i',coreness_CCS','type','Kendall');
% a1_core_burt=corr(average_i',coreness_burt','type','Kendall');
% %a11_core_burt=corr(average_i',coreness_newburt','type','Kendall');
% a2_neighbor_burt_core=corr(average_i',neighbor_burt_core','type','Kendall');
% a3_extend_neighbor_burt_core=corr(average_i',extend_neighbor_burt_core','type','Kendall');
% a4_neighbor_core=corr(average_i',neighbor_core','type','Kendall');
% a5_extend_neighbor_core=corr(average_i',extend_neighbor_core','type','Kendall');
