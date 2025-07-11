function cdt=core_discount_two(mixedsig,core,k)
[r,l]=size(mixedsig);
for i=1:r
 B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
 corelin=core;
 for j=1:length(B)
     C=find(mixedsig(B(j),:));
     D=intersect(B,C);%D表示两个节点之间的交集
%      if(core(i)>= core(B(j)))
     corelin(B(j))=core(B(j))*(1-length(D)/(k*k));
%      end
 end
 cdt(i)=sum(corelin(B));
end
    
end