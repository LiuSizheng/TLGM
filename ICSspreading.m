function [r_result,i_result]=ICSspreading(N,A,infection_node,irate,rrate)%A�����ڽӾ���N��������ڵ�������startnode��ʾ��ʼ��Ⱦ�ڵ㼯�ϣ�irate��Ⱦ����
    
%��ʼʱ�ڵ��״̬��,��ʼʱֻ�нڵ�1Ϊ��Ⱦ״̬�������Ķ�Ϊ�׸�Ⱦ״̬  
    start_node=infection_node;
    [m,n]=size(start_node);
    head=n;            
    %����β����ʼ����Ϊ�գ�tail==head
   tail=1;            
%��ͷ�м����ȾԴ�ڵ�
for i=1:n
queue(i)=start_node(i);      
%������չ
end

%��Ⱦ�ڵ��б� 
infection=start_node;  
%�ָ��ڵ��б�  
recover=[];
%�׸�Ⱦ�ڵ��б�


%�׸�Ⱦ�ڵ��б�
for i=1:N
    %��ʼʱ��start_nodeΪ��Ⱦ״̬
    if ismember(i,start_node)==1
        %-1��ʾ�ýڵ��Ѿ����б���ɾ��
        susceptible(i)=-1;
    else
    %��ʼʱ������start_nodeΪ��Ⱦ״̬�⣬�����ڵ㶼�����׸�Ⱦ״̬
     susceptible(i)=i;
    end
end

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
                 if infection_random < irate
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

