clear ;
 A=load('USAir.txt');
%       A=A+1;
% TT=A(:, 1:2);
% aas=TT(:,1)';
% bbs=TT(:,2)';
% GG=graph(aas,bbs);
% [edgebins,ac]=biconncomp(G);

% % %%%%%%%%%%%%%%%%%%数据处理%%%%%%%%%%%%%%%%%%%%%
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
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%数据处理%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TT=A(:, 1:2);
% mixedsig=zeros(max(max(TT)));
% len=length(TT);
% for i=1:len
%     mixedsig(TT(i,1),TT(i,2))=1;
%     mixedsig(TT(i,2),TT(i,1))=1;
% end

kk=sparse(mixedsig);
[a b]=components(kk);
[B] = largestcomponent(mixedsig);
mixedsig=mixedsig(B,B);

BV=mixedsig;
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

[a,b]=sort(core,'descend');
hebing=[a;b];

%下一步计算出每一个节点的邻居核数neighbor_core
m=mixedsig;
[r,l]=size(mixedsig);
neighbor_core=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_core(i)=0;
    for j=1:len1
        neighbor_core(i)=core(B(j))+neighbor_core(i);
    end
    
end
%接下来编写扩展算法extend_neighbor_core
[r,l]=size(mixedsig);
extend_neighbor_core=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    extend_neighbor_core(i)=0;
    for j=1:len1
        extend_neighbor_core(i)=neighbor_core(B(j))+extend_neighbor_core(i);
    end
end

C={};
[a,b]=sort(extend_neighbor_core,'descend');
[max_core,indexm]=max(extend_neighbor_core);%找出第i轮中得分最高的那个值
biggest_core=find(extend_neighbor_core==max_core);
len0=length(biggest_core);
DegreeAll=sum(mixedsig);
k=20;
K=200;  %%%节点的初始激活数目为K
for i=1:K
    aaa=biggest_core(find(DegreeAll(biggest_core)==max(DegreeAll(biggest_core))));
    len0=length(aaa);
    initial_set(i)=aaa(ceil(rand(1)*len0));%随机取得最大核数中的一个作为起始节点%%%%
    neighbor_of_i=find(mix(initial_set(i),:))%初始节点的邻居都找出来
    zz=find(core(neighbor_of_i)==core(initial_set(i)));%%找到初始节点的邻居中核数与其相等的节点
%     if isempty(zz)
%         zz=find(core(neighbor_of_i)==core(initial_set(i))-1);
%     end
    C{i}=neighbor_of_i(zz);%找到节点i的内部链接的各个节点,组成子图C
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
    [max_core,indexm]=max(core);%找出第i轮中得分最高的那个值
    biggest_core=find(core==max_core);
    len0=length(biggest_core);
    DegreeAll=sum(mix);
    
end
hebing_initial_number=[C_number_of_initial_i;initial_set]';
dfdf=sortrows(hebing_initial_number,'descend')';
new_initial_set=dfdf(2,1:k);


n=N;
%%%%%%%%%%%%%%%%%%k表示种子节点个数，n表示网络节点总数
%function state=VoteRank(mixedsig,k,n)
state=zeros(1,k);
DegreeAll=sum(mixedsig);
[B index]=sort(DegreeAll,'descend');%按度进行排序并记录下标
stack_degree=zeros(1,n);
stack_degree(index(1:k))=1;


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

stack_voterank=state;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%利用单纯的K-核排序方法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stack_k_shell=hebing(2,1:k);



bata=0.1;recover=1;run_number=1000;pp=1
for i=1:20
    disp(i);
    sa=state(1,1:pp);
    for m=1:run_number
        [r_result,i_result]=sir_mutiple_initial_set_weigthed_degree(N,BV,sa,recover);
        r1(m)=r_result;
    end
    result_VoteRank(i)=sum(r1)/run_number;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    sr=new_initial_set(1,1:pp);

    for m=1:run_number
        [r_result,i_result]=sir_mutiple_initial_set_weigthed_degree(N,BV,sr,recover);
        r2(m)=r_result;
    end
    result_run(i)=sum(r2)/run_number;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%利用单纯的度排序方法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sd=index(1,1:pp);
  
    for m=1:run_number
        [r_result,i_result]=sir_mutiple_initial_set_weigthed_degree(N,BV,sd,recover);
        r3(m)=r_result;
    end
    result_degree(i)=sum(r3)/run_number;

%     
    sk=hebing(2,1:pp);

    for m=1:run_number
        [r_result,i_result]=sir_mutiple_initial_set_weigthed_degree(N,BV,sk,recover);
        r4(m)=r_result;
    end
    result_k_shell(i)=sum(r4)/run_number;
 
    pp=pp+1;
end

result_all=[result_k_shell;result_degree;result_VoteRank;result_run];
bata_all=1:1:k;
figure;
axes('linewidth',1,'box','on','FontSize',16);
xlabel('感染概率','FontSize',16);
ylabel('感染比例','FontSize',16);
% grid on;
hold on;
plot(bata_all,result_k_shell,'o-',bata_all,result_degree,'v-',bata_all,result_VoteRank,'*-',bata_all,result_run,'p-');

legend('core','degree','VoteRank','run','Location','southeast');
% result_k_shell11=result_k_shell(1:5:k);
% result_degreellll=result_degree(1:5:k)
% result_VoteRankllll=result_VoteRank(1:5:k);
% result_runllll=result_run(1:5:k);
% result_all=[result_k_shell11;result_degreellll;result_VoteRankllll;result_runllll];
% bata_all=1:5:k;
% figure;
% axes('linewidth',1,'box','on','FontSize',16);
% xlabel('感染概率','FontSize',16);
% ylabel('感染比例','FontSize',16);
% % grid on;
% hold on;
% plot(bata_all,result_k_shell11,'o-',bata_all,result_degreellll,'v-',bata_all,result_VoteRankllll,'*-',bata_all,result_runllll,'p-');
% 
% 
% 
% % B={};k
% % for i=1:k
% %     a=i;
% %     b=C{i};
% 
% % 
% % bata=0.05;recover=1;run_number=100;
% % for i=1:6
% %     for m=1:run_number
% %         r1(m)=sir_more(A,stack_voterank,bata,recover);
% %     end
% %     result_VoteRank(i)=sum(r1)/run_number;
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     
% %     for m=1:run_number
% %         r2(m)=sir_more(A,stack_run,bata,recover);
% %     end
% %     result_run(i)=sum(r2)/run_number;
% %     
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%利用单纯的度排序方法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     
% %     for m=1:run_number
% %         r3(m)=sir_more(A,stack_degree,bata,recover);
% %     end
% %     result_degree(i)=sum(r3)/run_number;
% %     
% %     for m=1:run_number
% %         r4(m)=sir_more(A,stack_k_shell,bata,1);
% %     end
% %     result_k_shell(i)=sum(r4)/run_number;
% %     bata=bata+0.05;
% % end
% % 
