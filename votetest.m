clear ;
A=load('URVemail.txt');
%   A=A+1;
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
mixedsig=mixedsig(B,B);
A=mixedsig;
N=size(A,1);
b=sum(A,2);
k1=sum(b)/N;
k2=sum(b'*b)/N;
irate=k1/k2;
irate_real=0.15;
T=length(find(mixedsig==1));
T=T/2;
A=sparse(mixedsig);
core=core_numbers(A)';
[a,b]=sort(core,'descend');
hebing=[a;b];

n=N;k=50;
%%%%%%%%%%%%%%%%%%k表示种子节点个数，n表示网络节点总数
%function state=VoteRank(mixedsig,k,n)
state=zeros(1,k);
DegreeAll=sum(mixedsig);
[B index]=sort(DegreeAll,'descend');%按度进行排序并记录下标
VoteAbility=ones(1,n);%产生一个1行n列的数组，用于存放节点投票能力,最初每一个节点的投票能力都是1
state(1)=index(1);
VoteAbility(state(1))=0;
for i=2:k
    neigbor=find(mixedsig(state(i-1),:));
    VoteAbility(neigbor)= VoteAbility(neigbor)-1/k1;
    %%%将已经被选出来的节点的领域信息置0
    mixedsig(state(i-1),:)=0;
    mixedsig(:,state(i-1))=0;
    %%%计算出第i轮中所有节点的得分值
    for j=1:n
        score_allnode(j)=sum(VoteAbility(find(mixedsig(j,:))));
    end
    [maxa,indexm]=max(score_allnode);%找出第i轮中得分最高的那个值
    state(i)=indexm;
end
stack_voterank=zeros(1,n);
stack_voterank(state)=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%利用单纯的K-核排序方法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stack_k_shell=zeros(1,n);
stack_k_shell(hebing(2,1:k))=1;

%  for m=1:10000
%  r(m)=sir_more(A,stack_degree,bata,1);
%  end
%   result_degree=sum(r)/10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%利用单纯的度排序方法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stack_degree=zeros(1,n);
stack_degree(index(1:k))=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%利用单纯的K-核排序方法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  for m=1:10000
%  r(m)=sir_more(A,stack_degree,bata,1);
%  end
%   result_degree=sum(r)/10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=sparse(mixedsig);
core=core_numbers(A)';
[a,b]=sort(core,'descend');
hebing=[a;b];
stack_k_shell=zeros(1,n);
stack_k_shell(hebing(2,1:k))=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%我们的改进算法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VoteAbility_first=ones(1,n);%产生一个1行n列的数组，用于每一个节点对于直接邻居的投票能力,最初每一个节点的投票能力都是1
VoteAbility_first(1,:)=1;
%  VoteAbility_second=zeros(1,n);%产生一个1行n列的数组，用于每一个节点对于二级邻居节点的投票能力,最初每一个节点的投票能力都是1/(k1)^2
%  VoteAbility_second(1,:)=(1/k1)^2;
first_neighbor_set={};%节点的所有一阶邻居
second_neighbor_set={};%节点的所有二阶邻居
for p=1:n
    first_neighbor_set{p}=neighbor(A, p);
    s=length(first_neighbor_set{p});
    second_neighbor_set{p}=[];
    for i=1:s
        a=find(A(first_neighbor_set{p}(i),:));
        second_neighbor_set{p}=[second_neighbor_set{p},a];
    end
    second_neighbor_set{p}=unique(second_neighbor_set{p})
    second_neighbor_set{p}=setdiff(second_neighbor_set{p},[first_neighbor_set{p},p]);%找到节点p的二层邻居
end
for i=1:n
    statet(i)=length(first_neighbor_set{i})+length(second_neighbor_set{i})*(1/k1);
end
[maxb,indexb]=max(statet);
stack_proposed(1)=indexb;%第一个最大的找到了,将其双重投票能力置为0
VoteAbility_first(indexb)=0;
%  VoteAbility_second(indexb)=0;
for j=2:k
    
    VoteAbility_first(first_neighbor_set{stack_proposed(j-1)})=VoteAbility_first(first_neighbor_set{stack_proposed(j-1)})*(1/k1)^2;%%%将已经投过票的节点投票能力进行折扣
    VoteAbility_first(second_neighbor_set{stack_proposed(j-1)})=VoteAbility_first(second_neighbor_set{stack_proposed(j-1)})*(1/k1);%%%
    
    mixedsig(stack_proposed(j-1),:)=0;
    mixedsig(:,stack_proposed(j-1))=0;
    
    for i=1:n
        score_allnode(i)=sum(VoteAbility_first(first_neighbor_set{i}))+sum(VoteAbility_first(second_neighbor_set{i}))*(1/k1);
    end
    [maxa,indexm]=max(score_allnode);%找出第i轮中得分最高的那个值
    stack_proposed(j)=indexm;
end

ex_number=100;bata=0.05;

stack_voterank=zeros(1,n);
stack_voterank(state)=1;

stack_voterank_proposed=zeros(1,n);
stack_voterank_proposed(stack_proposed)=1;

stack_degree=zeros(1,n);
stack_degree(index(1:k))=1;
for i=1:6
    
for m=1:ex_number
    r1(m)=sir_more(A,stack_voterank,bata,1);
end
result_VoteRank(i)=sum(r1)/ex_number;

for m=1:ex_number
    r2(m)=sir_more(A,stack_degree,bata,1);
end
result_degree(i)=sum(r2)/ex_number;

for m=1:ex_number
    r3(m)=sir_more(A,stack_voterank_proposed,bata,1);
end
result_proposed(i)=sum(r3)/ex_number;

for m=1:ex_number
    r4(m)=sir_more(A,stack_k_shell,bata,1);
end
result_k_shell(i)=sum(r4)/ex_number;
bata=bata+0.05;
end


result_all=[result_k_shell;result_degree;result_VoteRank;result_proposed];
bata_all=[0.05,0.1,0.15,0.2,0.25,0.3];
figure;
axes('linewidth',1,'box','on','FontSize',16);
xlabel('感染概率','FontSize',16);
ylabel('感染比例','FontSize',16);
% grid on;
hold on;
plot(bata_all,result_k_shell,'o-',bata_all,result_degree,'v-',bata_all,result_VoteRank,'*-',bata_all,result_proposed,'p-');
legend('core','degree','VoteRank','RVoteRank');




















