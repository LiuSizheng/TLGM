
 function M_core=MDD(mixedsig)
b=sum(mixedsig);
bbbb=size(mixedsig,2);
tad=mixedsig;  %�ڽӾ���
len=length(mixedsig);  %���󳤶�
AA=[];
BB=[];
BB(1,2)=0;
    t=1;  %���Ƶ�j��  ����
 while(t==1)  %  ÿѭ��һ�Σ�tad�ı䣬ȥ������С��bbb�Ľڵ㣻ֱ��û�ж���С��bbb�Ľڵ�
     [SS,DD]=size(AA); 
     if DD==bbbb  %����Ԫ�ض������꣬�˳�
        break;
         end
        t=0;
        %Ѱ��ʣ�½ڵ����ۺ�ʣ�����С���Ǹ��㣬���ýڵ���ۺ϶���Ϊ��һ������Ŀ�ֵbbb
                 as=[];
                 for h=1:len
                     if sum(tad(h,:))==0&&(ismember(h,BB(:,2)))==1
                        
                         as(h)=1000;
                     else
                         as(h)=(b(h)-sum(tad(h,:)))*0.7+sum(tad(h,:));
                     end
                 end
             bbb=min(as);
        for i=1:len  %�Ӿ���1��len��, ����ÿһ���ڵ��ʣ��ȣ�����ڵ��Ѿ������룬��ô��ֱ��ָ���ýڵ�Ϊ��0
            if sum(tad(i,:))==0&&(ismember(i,BB(:,2)))==1
            k=sum(tad(i,:));%����i�ж���
            else
                k=(b(i)-sum(tad(i,:)))*0.7+sum(tad(i,:));
            end            
            if k==0  %����Ϊ0����һ��iֵ
                % t=1 ;  %
                continue;        
            elseif k<=bbb  %����С��j��
                t=1; %������һ�λ�Ҫѭ��
                tad(i,:)=0; %i�ڵ�ӵ�j�㣬��������Ϊ0������i��ֵ��Ϊ0
                tad(:,i)=0; %��Ӧi����Ϊ0
               AA=[AA;[bbb,i]];%��С�ڿ�ֵbbb�Ľڵ�������AA��        
               BB=AA;
                %for n=1:len   %�ж�������������Ϊ ��tad(:,i)Ԫ����Ϊ0 �������б�0
                  %  if sum(tad(n,:))==0&&(ismember(n,AA(:,2)))==0%���n�����������ģ���AA�е��Ҳ���һ�е�ֵ��Ӧ�ò�ͬ
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
    




