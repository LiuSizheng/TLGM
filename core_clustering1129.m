clear ;
A=load('infectious.txt');
%       A=A+1;
% TT=A(:, 1:2);
% aas=TT(:,1)';
% bbs=TT(:,2)';
% GG=graph(aas,bbs);
% [edgebins,ac]=biconncomp(G);

% % %%%%%%%%%%%%%%%%%%���ݴ���%%%%%%%%%%%%%%%%%%%%%
TT=A(:, 1:2);
b=size(TT,1);
hebing=[TT(:,1);TT(:,2)];
c=b*2;
hebing=hebing';


vv=unique(hebing);
cc=(1:length(vv));
duibi=[cc;vv];
for i=1:length(vv)
    r=find(hebing==duibi(2,i));
    l=length(r)
    for j=1:l
        hebing(r(j))=duibi(1,i);
    end
end
jieguo=[hebing(1:b);hebing((b+1):c)]';
TT=jieguo;
mixedsig=zeros(max(max(TT)));
len=length(TT);
for i=1:len
    mixedsig(TT(i,1),TT(i,2))=1;
    mixedsig(TT(i,2),TT(i,1))=1;
end
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%���ݴ���%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
mix=mixedsig;
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
cdist= kfindex(mixedsig);
jkf=sum(cdist);
cdist= kfindex(mixedsig);

for i=1:N
    B=find(mixedsig(i,:));%find�����������ؾ���mixedsig��i,:���з���Ԫ������λ��
    len1=length(B);%len1��ʾ����B�ĳ���,��ʾ�ڵ���ھ�����
    ddd2=[];
    for j=1:len1
        ddd2(j)=jkf(i)/cdist(i,B(j))*log(jkf(i)/cdist(i,B(j)));
    end
    gongshi8_2(i)=-1*sum(ddd2);
end


[a,b]=sort(core,'descend');
hebing=[a;b];

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

C={};
[a,b]=sort(extend_neighbor_core,'descend');
[max_core,indexm]=max(extend_neighbor_core);%�ҳ���i���е÷���ߵ��Ǹ�ֵ
biggest_core=find(extend_neighbor_core==max_core);
len0=length(biggest_core);
DegreeAll=sum(mixedsig);
k=20;
K=20;  %%%�ڵ�ĳ�ʼ������ĿΪK
for i=1:K
    aaa=biggest_core(find(DegreeAll(biggest_core)==max(DegreeAll(biggest_core))));
    len0=length(aaa);
    initial_set(i)=aaa(ceil(rand(1)*len0));%���ȡ���������е�һ����Ϊ��ʼ�ڵ�%%%%
    neighbor_of_i=find(mix(initial_set(i),:))%��ʼ�ڵ���ھӶ��ҳ���
    zz=find(core(neighbor_of_i)>=core(initial_set(i)));
    C{i}=neighbor_of_i(zz);%�ҵ��ڵ�i���ڲ����ӵĸ����ڵ�,�����ͼC
    C{i}=[C{i},initial_set(i)];
    first_neighbor_of_C=neighbor_of_node_set(mix,C{i});
    degree_of_first_neighbor_of_C=DegreeAll(first_neighbor_of_C);
    togther=[degree_of_first_neighbor_of_C;first_neighbor_of_C]';
    togtherr=sortrows(togther)';
    for j=1:length(first_neighbor_of_C)
        ddd=neighbor(mix,togtherr(2,j));
        jiaoji=intersect(ddd,C{i});
        chaji=setdiff(ddd,jiaoji);
        if length(jiaoji)>length(chaji)
            C{i}=[C{i},togtherr(2,j)];
        end
        
    end
    
    second_neighbor_of_C=neighbor_of_node_set(mix,C{i});
    degree_of_second_neighbor_of_C=DegreeAll(second_neighbor_of_C);
    togther=[degree_of_second_neighbor_of_C;second_neighbor_of_C]';
    togtherr=sortrows(togther)';
    for j=1:length(second_neighbor_of_C)
        ddd=neighbor(mix,togtherr(2,j));
        jiaoji=intersect(ddd,C{i});
        chaji=setdiff(ddd,jiaoji);
        if length(jiaoji)>length(chaji)
            C{i}=[C{i},togtherr(2,j)];
        end
        
    end
    third_neighbor_of_C=neighbor_of_node_set(mix,C{i});
    degree_of_third_neighbor_of_C=DegreeAll(third_neighbor_of_C);
    togther=[degree_of_third_neighbor_of_C;third_neighbor_of_C]';
    togtherr=sortrows(togther)';
    for j=1:length(third_neighbor_of_C)
        ddd=neighbor(mix,togtherr(2,j));
        jiaoji=intersect(ddd,C{i});
        chaji=setdiff(ddd,jiaoji);
        if length(jiaoji)>length(chaji)
            C{i}=[C{i},togtherr(2,j)];
        end
        
    end
    C_number_of_initial_i(i)=length(C{i});
    linshi3=C{i};
    
    linshi4=linshi3(find(DegreeAll(linshi3)==max(DegreeAll(linshi3))));
    initial_set(i)=linshi4(ceil(rand(1)*length(linshi4)));
    
    mix(C{i},:)=0;
    mix(:,C{i})=0;
    [a,b]=sort(core,'descend');
    MN=sparse(mix);
    core=core_numbers(MN)';
    [a,b]=sort(core,'descend');
    [max_core,indexm]=max(core);%�ҳ���i���е÷���ߵ��Ǹ�ֵ
    biggest_core=find(core==max_core);
    len0=length(biggest_core);
    DegreeAll=sum(mix);
    
end
hebing_initial_number=[C_number_of_initial_i;initial_set]';
dfdf=sortrows(hebing_initial_number,-1)';
new_initial_set=dfdf(2,1:k);
stack_run=zeros(1,N);
stack_run(new_initial_set)=1;




%%%%%%%%%%%%%%%%%%

n=N;
%%%%%%%%%%%%%%%%%%k��ʾ���ӽڵ������n��ʾ����ڵ�����
%function state=VoteRank(mixedsig,k,n)
state=zeros(1,k);
DegreeAll=sum(mixedsig);
[B index]=sort(DegreeAll,'descend');%���Ƚ������򲢼�¼�±�
stack_degree=zeros(1,n);
stack_degree(index(1:k))=1;


VoteAbility=ones(1,n);%����һ��1��n�е����飬���ڴ�Žڵ�ͶƱ����,���ÿһ���ڵ��ͶƱ��������1
state(1)=index(1);
VoteAbility(state(1))=0;
for i=2:k
    neigbor=find(mixedsig(state(i-1),:));
    VoteAbility(neigbor)= VoteAbility(neigbor)-1/k1;
    %%%���Ѿ���ѡ�����Ľڵ��������Ϣ��0
    mixedsig(state(i-1),:)=0;
    mixedsig(:,state(i-1))=0;
    %%%�������i�������нڵ�ĵ÷�ֵ
    for j=1:n
        score_allnode(j)=sum(VoteAbility(find(mixedsig(j,:))));
    end
    [maxa,indexm]=max(score_allnode);%�ҳ���i���е÷���ߵ��Ǹ�ֵ
    state(i)=indexm;
end
stack_voterank=zeros(1,n);
stack_voterank(state)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���õ�����K-�����򷽷�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stack_k_shell=zeros(1,n);
stack_k_shell(hebing(2,1:k))=1;
   [B index]=sort(gongshi8_2,'ascend');%���Ƚ������򲢼�¼�±�
stack_8_2=zeros(1,n);
stack_8_2(index(1:k))=1;

bata=0.05;recover=1;run_number=100;
for i=1:6
    for m=1:run_number
        r1(m)=sir_more(A,stack_voterank,bata,recover);
    end
    result_VoteRank(i)=sum(r1)/run_number;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for m=1:run_number
        r2(m)=sir_more(A,stack_run,bata,recover);
    end
    result_run(i)=sum(r2)/run_number;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���õ����Ķ����򷽷�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for m=1:run_number
        r3(m)=sir_more(A,stack_degree,bata,recover);
    end
    result_degree(i)=sum(r3)/run_number;
    
    for m=1:run_number
        r4(m)=sir_more(A,stack_k_shell,bata,1);
    end
    result_k_shell(i)=sum(r4)/run_number;
    
    for m=1:run_number
        r5(m)=sir_more(A,stack_8_2,bata,1);
    end
    result_8_2(i)=sum(r5)/run_number;
    bata=bata+0.05;
end

result_all=[result_k_shell;result_degree;result_VoteRank;result_run;result_8_2];
bata_all=[0.05,0.1,0.15,0.2,0.25,0.3];
figure;
axes('linewidth',1,'box','on','FontSize',16);
xlabel('��Ⱦ����','FontSize',16);
ylabel('��Ⱦ����','FontSize',16);
% grid on;
hold on;
plot(bata_all,result_k_shell,'o-',bata_all,result_degree,'v-',bata_all,result_VoteRank,'*-',bata_all,result_run,'p-',bata_all,result_8_2,'*-');
legend('core','degree','VoteRank','RVoteRank','Location','southeast','result_8_2');



ll=[];B={};new=[];
for i=1:k
    a=i;
    b=C{i};
    aa=repmat(a,1,length(C{i}));
   ll=[b',aa'];
   B{i}=ll;
   new=[new;ll];
end

%    stack_ac=zeros(1,n);
%  stack_ac(ic(1:k))=1;
%
%  for m=1:100
%  r(m)=sir_more(A,stack_ac,bata,recover);
%  end
%   result_ac=sum(r)/100;
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%