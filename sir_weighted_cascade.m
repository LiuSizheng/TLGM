function r=sir_weighted_cascade(w,state,mu)
%�������w����ʼ״̬state����Ⱦ��bata���ظ���mu
%�����Ⱦ����r

n=length(w);
%0 S̬��1 I̬��2 R̬



%��Ⱦ���� һ���ڵ�һ���ڵ��Ⱦ
while(sum(state==1))            %�ڵ�û��1̬����
    state_i=find(state==1);     %Ѱ�ҽڵ�Ϊ1�Ľڵ�
    nn=length(state_i);         %1̬�ĸ���
%     rand_2=rand(2,nn);          %���������        
    for i=1:nn
        if rand(1)<mu      %I̬���R̬
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