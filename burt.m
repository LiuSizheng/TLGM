%����ÿһ���ڵ��burtϵ��
function cc=burt(mixedsig) 

  [r,l]=size(mixedsig);
  for i=1:r
       B=find(mixedsig(i,:));%find�����������ؾ���mixedsig��i,:���з���Ԫ������λ��
       len1=length(B);%len1��ʾ����B�ĳ���,��ʾ�ڵ���ھ�����
       com_neighbor=cell(len1,1);%com_neighbor�����洢�ڵ�i��ÿһ���ھ�
       for j=1:len1
          
           nei=find(mixedsig(B(j),:));
           num=length(nei);%num��ʾ�ڵ��ھӵĶ�
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