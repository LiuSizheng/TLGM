function cdt=core_discount(mixedsig,core,k)
[r,l]=size(mixedsig);
for i=1:r
 B=find(mixedsig(i,:));%find�����������ؾ���mixedsig��i,:���з���Ԫ������λ��
 
 for j=1:length(B)
     C=find(mixedsig(B(j),:));
     D=intersect(B,C);%D��ʾ�����ڵ�֮��Ľ���
     corelin=core;
     corelin(B(j))=core(B(j))*(1-length(D)/(k*k));
 end
 cdt(i)=sum(corelin(B));
end
    
end