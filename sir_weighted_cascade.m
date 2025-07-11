function r=sir_weighted_cascade(w,state,mu)
%输入矩阵w；初始状态state；传染率bata；回复率mu
%输出感染比例r

n=length(w);
%0 S态，1 I态，2 R态



%感染过程 一个节点一个节点感染
while(sum(state==1))            %节点没有1态结束
    state_i=find(state==1);     %寻找节点为1的节点
    nn=length(state_i);         %1态的个数
%     rand_2=rand(2,nn);          %产生随机数        
    for i=1:nn
        if rand(1)<mu      %I态变成R态
             state(state_i(i))=2;
        end
          lj=find(w(state_i(i),:)==1);
          for j=1:length(lj)
              if state(lj(j))==0&&rand(1)<(1/degree(lj(j),w))
                  state(lj(j))=1;
              end
          end
        
    end
end
r=length(find(state==2));