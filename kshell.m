function ret_r=kshell(mixedsig)
j=1;   %��j��
t=1;  %
r=1;
tad=mixedsig;  %�ڽӾ���
k=0;
ii=1;
ret_m=[];  %����ֵ��  i��j  �� ��i�㣺���нڵ�
ret_r=[]; 
len=length(mixedsig);  %���󳤶�
tt=1;
n=1;
nn=1;
m_t=[];
m_tt=[];
m_ttt=[];
p=1;
while (tt==1)   %������߲����
  
    if sum(sum(tad))==0  %����Ԫ��Ϊ0�����˳�
        break;
    end
    t=1;  %���Ƶ�j��  ����
    while(t==1)  %  ÿѭ��һ�Σ�tad�ı䣬ȥ������С��j�Ľڵ㣻ֱ��û�ж���С��j�Ľڵ�
        t=0;
        ii=1; %��j���ii���ڵ�
        for i=1:len  %�Ӿ���1��len��,  ȥ��С��j�Ľڵ�
            k=sum(tad(i,:));   %����i�ж���
            if k==0  %����Ϊ0����һ��iֵ
                % t=1 ;  %
                continue;
            elseif k<=j  %����С��j��
                t=1; %������һ�λ�Ҫѭ��
                tad(i,:)=0;  %i�ڵ�ӵ�j�㣬��������Ϊ0������i��ֵ��Ϊ0
                tad(:,i)=0; %��Ӧi����Ϊ0

                ret_m(j,ii)=i;   %��i�ڵ�ӵ�j��
               ret_r(j,r)=i;
               r=r+1;
                m_ttt=union(ret_m(j,ii),m_ttt);
                for n=1:len   %�ж�������������Ϊ ��tad(:,i)Ԫ����Ϊ0 �������б�0
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