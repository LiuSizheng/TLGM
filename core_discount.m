function cdt=core_discount(mixedsig,core,k)
[r,l]=size(mixedsig);
for i=1:r
 B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
 
 for j=1:length(B)
     C=find(mixedsig(B(j),:));
     D=intersect(B,C);%D表示两个节点之间的交集
     corelin=core;
     corelin(B(j))=core(B(j))*(1-length(D)/(k*k));
 end
 cdt(i)=sum(corelin(B));
end
    
end