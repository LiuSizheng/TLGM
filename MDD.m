
 function M_core=MDD(mixedsig)
b=sum(mixedsig);
bbbb=size(mixedsig,2);
tad=mixedsig;  %邻接矩阵
len=length(mixedsig);  %矩阵长度
AA=[];
BB=[];
BB(1,2)=0;
    t=1;  %控制第j层  计算
 while(t==1)  %  每循环一次，tad改变，去掉度数小于bbb的节点；直到没有度数小于bbb的节点
     [SS,DD]=size(AA); 
     if DD==bbbb  %所有元素都排序完，退出
        break;
         end
        t=0;
        %寻找剩下节点中综合剩余度最小的那个点，将该节点的综合度作为下一步剥离的壳值bbb
                 as=[];
                 for h=1:len
                     if sum(tad(h,:))==0&&(ismember(h,BB(:,2)))==1
                        
                         as(h)=1000;
                     else
                         as(h)=(b(h)-sum(tad(h,:)))*0.7+sum(tad(h,:));
                     end
                 end
             bbb=min(as);
        for i=1:len  %从矩阵1至len行, 计算每一个节点的剩余度，如果节点已经被剥离，那么就直接指定该节点为度0
            if sum(tad(i,:))==0&&(ismember(i,BB(:,2)))==1
            k=sum(tad(i,:));%计算i行度数
            else
                k=(b(i)-sum(tad(i,:)))*0.7+sum(tad(i,:));
            end            
            if k==0  %度数为0，下一个i值
                % t=1 ;  %
                continue;        
            elseif k<=bbb  %度数小于j层
                t=1; %控制下一次还要循环
                tad(i,:)=0; %i节点加到j层，将度数至为0，所有i行值至为0
                tad(:,i)=0; %相应i列至为0
               AA=[AA;[bbb,i]];%将小于壳值bbb的节点剥离放在AA中        
               BB=AA;
                %for n=1:len   %判断其他行有无因为 将tad(:,i)元素设为0 而所有行变0
                  %  if sum(tad(n,:))==0&&(ismember(n,AA(:,2)))==0%这个n必须是新增的，和AA中的右侧那一列的值都应该不同
                  %     if (b(n)-sum(tad(n,:)))*0.7+sum(tad(h,:));
                 %      AA=[AA;[bbb,n]];
                  %     AA=unique(AA,'rows');
                  %  end
                end %for n=1:len
        end % if k==0
end %end of  i=1:len
         
 AA=sortrows(AA,2)   ;
    AA=AA';
M_core=AA(1,:);
 end
    




