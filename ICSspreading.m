function [r_result,i_result]=ICSspreading(N,A,infection_node,irate,rrate)%A代表邻接矩阵，N代表网络节点总数，startnode表示起始感染节点集合，irate感染概率
    
%初始时节点的状态表,初始时只有节点1为感染状态，其他的都为易感染状态  
    start_node=infection_node;
    [m,n]=size(start_node);
    head=n;            
    %队列尾，开始队列为空，tail==head
   tail=1;            
%向头中加入感染源节点
for i=1:n
queue(i)=start_node(i);      
%队列扩展
end

%感染节点列表 
infection=start_node;  
%恢复节点列表  
recover=[];
%易感染节点列表


%易感染节点列表
for i=1:N
    %初始时，start_node为感染状态
    if ismember(i,start_node)==1
        %-1表示该节点已经从列表中删除
        susceptible(i)=-1;
    else
    %初始时，除了start_node为感染状态外，其他节点都处于易感染状态
     susceptible(i)=i;
    end
end

%开始按照广度优先搜索顺序向邻居节点传播
%判断队列是否为空
while tail~=head   
    %取队尾节点 
    i=queue(tail);  
    %如果该节点不在移除列表之中
    if isempty(find(recover==i,1))
       for j=1:N
             %如果节点j与当前节点i相连并且节点j不在感染列表与恢复列表中
            if (A(i,j)==1 && isempty(find(infection==j,1)) && isempty(find(recover==j,1)))
                 infection_random=rand(1);
                 if infection_random < irate
                    %新节点入列
                    queue(head)=j;  
                    %扩展队列
                    head=head+1;   
                    %将新节点j加入感染列表
                    infection=[infection j]; 
                    
                    %从易感染节点列表中删除该节点,设置为-1
                    [~,col,v] = find(susceptible==j) ;
                    susceptible(col)=-1;
                    susceptible(susceptible==-1)=[];                    
                 end
            end
        end
        %将感染的节点按概率加入恢复节点列表  
        recover_random=rand(1);
        if recover_random < rrate
            %恢复
            recover=[recover i];  
            %从感染列表中删除
            [~,col,v] = find(infection==i) ;
            infection(col)=-1;
            infection(infection==-1)=[];
             tail=tail+1; 
        end
       
        
    end %end if  isempty(find(recover==i,1)
end %end while
i_result=length(infection);
r_result=length(recover);

end

