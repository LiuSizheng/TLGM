function [r_result,i_result]=sir_mutiple_initial_set_weigthed_degree(N,A,infection_node,rrate)%A�����ڽӾ���N��������ڵ�������startnode��ʾ��ʼ��Ⱦ�ڵ㣬irate��Ⱦ����
    
%��ʼʱ�ڵ��״̬��,��ʼʱֻ�нڵ�1Ϊ��Ⱦ״̬�������Ķ�Ϊ�׸�Ⱦ״̬  
len=length(infection_node);
    start_node=infection_node(len);
         
    %����β����ʼ����Ϊ�գ�tail==head
   tail=1;            
%��ͷ�м����ȾԴ�ڵ�
queue(1:len)=infection_node;

%������չ
head=len+1;  

%��Ⱦ�ڵ��б� 
infection=infection_node;  
%�ָ��ڵ��б�  
recover=[];
%�׸�Ⱦ�ڵ��б�
susceptible=1:1:N;
susceptible(infection)=-1;


%��ʼ���չ����������˳�����ھӽڵ㴫��
%�ж϶����Ƿ�Ϊ��
while tail~=head   
    %ȡ��β�ڵ� 
    i=queue(tail);  
    %����ýڵ㲻���Ƴ��б�֮��
    if isempty(find(recover==i,1))
       for j=1:N
             %����ڵ�j�뵱ǰ�ڵ�i�������ҽڵ�j���ڸ�Ⱦ�б���ָ��б���
            if (A(i,j)==1 && isempty(find(infection==j,1)) && isempty(find(recover==j,1)))
                 infection_random=rand(1);
                 if infection_random < 1/degree(j,A)
                    %�½ڵ�����
                    queue(head)=j;  
                    %��չ����
                    head=head+1;   
                    %���½ڵ�j�����Ⱦ�б�
                    infection=[infection j]; 
                    
                    %���׸�Ⱦ�ڵ��б���ɾ���ýڵ�,����Ϊ-1
                    [~,col,v] = find(susceptible==j) ;
                    susceptible(col)=-1;
                    susceptible(susceptible==-1)=[];                    
                 end
            end
        end
        %����Ⱦ�Ľڵ㰴���ʼ���ָ��ڵ��б�  
        recover_random=rand(1);
        if recover_random < rrate
            %�ָ�
            recover=[recover i];  
            %�Ӹ�Ⱦ�б���ɾ��
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

