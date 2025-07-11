%计算每一个节点的burt系数
function cc=burt(mixedsig) 

  [r,l]=size(mixedsig);
  for i=1:r
       B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
       len1=length(B);%len1表示数组B的长度,表示节点的邻居总数
       com_neighbor=cell(len1,1);%com_neighbor用来存储节点i的每一个邻居
       for j=1:len1
          
           nei=find(mixedsig(B(j),:));
           num=length(nei);%num表示节点邻居的度
           com_neighbor{j}=intersect(find(mixedsig(i,:)),find(mixedsig(B(j),:)));
 totol_nei=length(com_neighbor{j});
 s=[];
 for k=1:totol_nei
     s(k)=0;
     a=1/degree(com_neighbor{j}(k),mixedsig);
     b=1/degree(i,mixedsig);
     s(k)=(a)*(b)+s(k);
 end
 cv(j)=(1/degree(i,mixedsig)+sum(s))^2;
       end
       cc(i)=sum(cv);
       cv=[];
  end
end