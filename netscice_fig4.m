clear ;
A=load('netscience.txt');
%A=A+1;
r=0.2;
TT=A(:, 1:2);
mixedsig=zeros(max(max(TT)));
N=size(mixedsig,1);
len=length(TT);
for i=1:len
    mixedsig(TT(i,1),TT(i,2))=1;
    mixedsig(TT(i,2),TT(i,1))=1;
end
kk=sparse(mixedsig);
[a b]=components(kk);
[B]=largestcomponent(mixedsig);
mixedsig=mixedsig(B,B);
A=mixedsig;
N=size(A,1);
b=sum(A,2);
k1=sum(b)/N;
k2=sum(b'*b)/N;
irate=k1/k2;
hefei1={};
for i=1:20%%%%%%r从0.05到1
    disp(i);
    hefei1{i}= secondfig(mixedsig,0.05*i);
end
load('netscience_average_sir1_26.mat', 'average_iall') ;
ss=[];
for j=1:26%%%%%%传播率从0.01到0.2
    for i=1:20%%%%%%r从0.05到1
        a00_hefei(i)=corr(average_iall{j}',hefei1{i}','type','Kendall');
    end
    ss(j,:)=a00_hefei;
end
x=0.06:0.01:0.20;
y=0.05:0.05:1;
%yy=fliplr(y); 
pd = makedist('Normal');
value=ss(6:20,:);
%value=ss;
tyy=value';
image(x,y,value','CDataMapping','scaled');%后面的参数是根据数值大小选取colorbar着色范围

set(gca,'YDir','normal'); %若y轴为时间，可将y轴倒置

colormap jet;colorbar;


ft=ss';
xlabel('β');ylabel('λ');

