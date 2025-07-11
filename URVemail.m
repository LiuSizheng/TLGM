 
clear all;
%irate=0.15;
%A=load('erdos.txt');
A=load('URVemail.txt');
%A=load('blog.txt');


 TT=A(:, 1:2);
 mixedsig=zeros(max(max(TT)));
 len=length(TT);
 for i=1:len
     mixedsig(TT(i,1),TT(i,2))=1;
     mixedsig(TT(i,2),TT(i,1))=1;
 end
kk=sparse(mixedsig);
[a b]=components(kk);
 [B] = largestcomponent(mixedsig);
%kk=sparse(mixedsig);
%[a b]=components(kk);
% [B] = largestcomponent(mixedsig);
mixedsig=mixedsig(B,B);
 A=mixedsig;
  N=size(A,1); 
  b=sum(A,2);
  k1=sum(b)/N;
  k2=sum(b'*b)/N; 
irate=k1/k2;
average_i=sir_main(mixedsig,500) ;
 %����������ڵ�ļ���ϵ��
 N=size(mixedsig,2);
C=zeros(1,N);
for i=1:N
    a=find(mixedsig(i,:)==1); %Ѱ����ͼ���ھӽڵ�
    b=find(mixedsig(:,i)==1);
    m=union(a,b'); 
    k=length(m);
    if k==1
        disp(['�ڵ�',int2str(i),'ֻ��һ���ھӽڵ㣬�����ϵ��Ϊ0']);
        C(i)=0;
    else
        B=mixedsig(m,m);
        C(i)=length(find(B==1))/(k*(k-1));
    end
end
aver_C=mean(C);
max_C=max(C);


%average_i=[5.29000000000000	4.88000000000000	4.780000000000000	4.53000000000000	3.90000000000000	3.79000000000000	3.05000000000000	3.04000000000000	2.50000000000000	2.51000000000000	2.09000000000000	2.05000000000000	2.02000000000000	2.34000000000000	1.61000000000000	4.34000000000000	2.19000000000000	2.22000000000000	2.21000000000000	2.24000000000000	2.23000000000000	2.20000000000000	2.67000000000000	2.32000000000000	1.56000000000000	1.55000000000000
%];
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
    disp(tt);
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
%disp(ret_m);


%����������������������core��
core=[];
[r,l]=size(ret_r);
for i=1:r
B=find(ret_r(i,:));
k1=length(B);
for j=1:k1
    core(ret_r(i,B(j)))=i;
end
end

 
 
 %�Ľ��㷨coreness_burt
 N=size(mixedsig,2);
 cc=burt(mixedsig);
 for i=1:N
     %extend_cluster_coreness(i)=core(i)+neighbor_core(i)/power(1.5,C(i));
     coreness_burt(i)=core(i)/power(10,cc(i));
 end
 %��дcoreness-burt����չ�㷨
 [r,l]=size(mixedsig);
extend_coreness_burt=[];
 for i=1:r
 B=find(mixedsig(i,:));%find�����������ؾ���mixedsig��i,:���з���Ԫ������λ��
    len1=length(B);%len1��ʾ����B�ĳ���
   extend_coreness_burt(i)=0;
     for j=1:len1  
    extend_coreness_burt(i)=coreness_burt(B(j))+extend_coreness_burt(i);
     end
 end
 %��һ�������ÿһ���ڵ���ھӺ���neighbor_core
m=mixedsig;
[r,l]=size(mixedsig);
neighbor_core=[];
for i=1:r
    B=find(mixedsig(i,:));%find�����������ؾ���mixedsig��i,:���з���Ԫ������λ��
    len1=length(B);%len1��ʾ����B�ĳ���
    neighbor_core(i)=0;
     for j=1:len1  
    neighbor_core(i)=core(B(j))+neighbor_core(i);
    end
    
end
%��������д��չ�㷨extend_neighbor_core
 [r,l]=size(mixedsig);
extend_neighbor_core=[];
 for i=1:r
 B=find(mixedsig(i,:));%find�����������ؾ���mixedsig��i,:���з���Ԫ������λ��
    len1=length(B);%len1��ʾ����B�ĳ���
   extend_neighbor_core(i)=0;
     for j=1:len1  
    extend_neighbor_core(i)=neighbor_core(B(j))+extend_neighbor_core(i);
     end
 end 


%��һ��������Ľ����ÿһ���ڵ���ھӺ���neighbor_burt_core
m=mixedsig;
[r,l]=size(mixedsig);
neighbor_burt_core=[];
for i=1:r
    B=find(mixedsig(i,:));%find�����������ؾ���mixedsig��i,:���з���Ԫ������λ��
    len1=length(B);%len1��ʾ����B�ĳ���
    neighbor_burt_core(i)=0;
     for j=1:len1  
    neighbor_burt_core(i)=coreness_burt(B(j))+neighbor_burt_core(i);
    end
    
end
%�Ľ������չ�㷨extend_neighbor_burt_core
 [r,l]=size(mixedsig);
extend_neighbor_burt_core=[];
 for i=1:r
 B=find(mixedsig(i,:));%find�����������ؾ���mixedsig��i,:���з���Ԫ������λ��
    len1=length(B);%len1��ʾ����B�ĳ���
   extend_neighbor_burt_core(i)=0;
     for j=1:len1  
    extend_neighbor_burt_core(i)=neighbor_burt_core(B(j))+extend_neighbor_burt_core(i);
     end
 end
qqqq0_core=corr(average_i',core','type','Kendall');
qqqq1_core_burt=corr(average_i',coreness_burt','type','Kendall');
qqqq2_neighbor_burt_core=corr(average_i',neighbor_burt_core','type','Kendall');
qqqq3_extend_neighbor_burt_core=corr(average_i',extend_neighbor_burt_core','type','Kendall');
qqqq4_neighbor_core=corr(average_i',neighbor_core','type','Kendall');
qqqq5_extend_neighbor_core=corr(average_i',extend_neighbor_core','type','Kendall');
 