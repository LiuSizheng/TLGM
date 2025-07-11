function ret_r=kshell(mixedsig)
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
end
%disp(ret_m);