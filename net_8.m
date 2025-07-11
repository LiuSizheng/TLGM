clear ;
A=load('rt-twitter-copen.txt');
    A=A+1;
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
% cc_ij=out_edge(mixedsig,3)
 A=mixedsig;
%  A=polblogs;
%  mixedsig=polblogs;
  N=size(A,1); 
  b=sum(A,2);
  k1=sum(b)/N;
  k2=sum(b'*b)/N; 
irate=k1/(k2-0);
core=core_numbers(sparse(mixedsig))';
bet=betweenness_centrality(sparse(mixedsig));
clo=Closeness_Centrality(mixedsig);

N=size(mixedsig,1);
d = degreea( mixedsig );kre=sum(d)/N;
% cc_ij=out_edge(mixedsig,3)
A=mixedsig;
  [D,aver_D]=Aver_Path_Length(mixedsig);
d = degreea( mixedsig );
core=core_numbers(sparse(mixedsig))';
maxcoree=max(core);
adddd=sum(d)/N;
CCER=sum(clustering_coefficients(sparse(mixedsig)))/N;
maxd=max(max(D));
acoree=sum(core)/N;
coef=adddd/acoree;


for i=1:N
    disp(i);
    j=1;%第一层
    nei0=[i];%第0层
    nei1=find(mixedsig(i,:));%算出节点i的直接邻居
    par=[];
    nei1=find(mixedsig(i,:));
    len1=length(nei1);
    par1=[];
    for j=1:len1
        par1(j)=d(i)*d(nei1(j));
    end
    
    par11=sum(par1);
    j=2;%第二层
    %     nei2=neighbor(mixedsig,nei1');
    % a=neighbor(mixedsig,nei1);
    a=[];
    for v=1:length(nei1)
        %     rr=mixedsig(nei1(v));
        k=find(mixedsig(nei1(v),:));
        a=union(k,a);
    end
    c=setdiff(a,nei0);%得到节点i的二跳邻居
    nei2=setdiff(c,nei1);%得到节点i的二跳邻居
    
    par2=[];
    for k=1:length(nei2)
        par2(k)=d(i)*d(nei2(k))/2^2;
        %         par2(k)=BBB(i,nei2(k))*(1/k1)^2+BBA(i,nei2(k))*(1/k1)^3;
    end
    par22=sum(par2);
    
    j=3;%第3层
    a=[];
    for v=1:length(nei2)
        %     rr=mixedsig(nei1(v));
        k=find(mixedsig(nei2(v),:));
        a=union(k,a);
    end
    c=setdiff(a,union(nei0,nei1));%得到节点i的3跳邻居
    nei3=setdiff(c,nei2);
    % BBA4=mixedsig^4;
    par3=[];
    for k=1:length(nei3)
        par3(k)=d(i)*d(nei3(k))/3^2;
    end
    par33=sum(par3);
    rund(i)=par11+par22+par33;
%      rund(i)=par11;
end


cij=max(core)-min(core);
for i=1:N
    disp(i);
    j=1;%第一层
    nei0=[i];%第0层
    nei1=find(mixedsig(i,:));%算出节点i的直接邻居
    par=[];
    nei1=find(mixedsig(i,:));
    len1=length(nei1);
    par1=[];
    for j=1:len1
        par1(j)=exp((core(i)-core(nei1(j)))/cij)*d(i)*d(nei1(j));
    end
    
    par11=sum(par1);
    j=2;%第二层
    %     nei2=neighbor(mixedsig,nei1');
    % a=neighbor(mixedsig,nei1);
    a=[];
    for v=1:length(nei1)
        %     rr=mixedsig(nei1(v));
        k=find(mixedsig(nei1(v),:));
        a=union(k,a);
    end
    c=setdiff(a,nei0);%得到节点i的二跳邻居
    nei2=setdiff(c,nei1);%得到节点i的二跳邻居
    
    par2=[];
    for k=1:length(nei2)
        par2(k)=exp((core(i)-core(nei2(k)))/cij)*d(i)*d(nei2(k))/2^2;
        %         par2(k)=BBB(i,nei2(k))*(1/k1)^2+BBA(i,nei2(k))*(1/k1)^3;
    end
    par22=sum(par2);
    
    j=3;%第3层
    a=[];
    for v=1:length(nei2)
        %     rr=mixedsig(nei1(v));
        k=find(mixedsig(nei2(v),:));
        a=union(k,a);
    end
    c=setdiff(a,union(nei0,nei1));%得到节点i的3跳邻居
    nei3=setdiff(c,nei2);
    % BBA4=mixedsig^4;
    par3=[];
    for k=1:length(nei3)
        par3(k)=exp((core(i)-core(nei3(k)))/cij)*d(i)*d(nei3(k))/3^2;
    end
    par33=sum(par3);
    runim(i)=par11+par22+par33;
%      runim(i)=par11;
end

ccg=burt(mixedsig);
coefda=acoree/adddd;



[r,l]=size(mixedsig);

for i=1:N
    disp(i);
    j=1;%第一层
    nei0=[i];%第0层
    nei1=find(mixedsig(i,:));%算出节点i的直接邻居
    par=[];
    nei1=find(mixedsig(i,:));
    len1=length(nei1);
    par1=[];
    for j=1:len1
        par1(j)=exp(-1*ccg(i))*(coefda*d(i)+core(i))*(coefda*d(nei1(j))+core(nei1(j)));
    end
    
    par11=sum(par1);
    j=2;%第二层
    %     nei2=neighbor(mixedsig,nei1');
    % a=neighbor(mixedsig,nei1);
    a=[];
    for v=1:length(nei1)
        %     rr=mixedsig(nei1(v));
        k=find(mixedsig(nei1(v),:));
        a=union(k,a);
    end
    c=setdiff(a,nei0);%得到节点i的二跳邻居
    nei2=setdiff(c,nei1);%得到节点i的二跳邻居
    
    par2=[];
    for k=1:length(nei2)
        par2(k)=exp(-1*ccg(i))*(coefda*d(i)+core(i))*(coefda*d(nei2(k))+core(nei2(k)))/2^2;
        %         par2(k)=BBB(i,nei2(k))*(1/k1)^2+BBA(i,nei2(k))*(1/k1)^3;
    end
    par22=sum(par2);
    
    j=3;%第3层
    a=[];
    for v=1:length(nei2)
        %     rr=mixedsig(nei1(v));
        k=find(mixedsig(nei2(v),:));
        a=union(k,a);
    end
    c=setdiff(a,union(nei0,nei1));%得到节点i的3跳邻居
    nei3=setdiff(c,nei2);
    % BBA4=mixedsig^4;
    par3=[];
    for k=1:length(nei3)
        par3(k)=exp(-1*ccg(i))*(coefda*d(i)+core(i))*(coefda*d(nei3(k))+core(nei3(k)))/3^2;
    end
    par33=sum(par3);
    run1(i)=par11+par22+par33;
%      run1(i)=par11;
end

extend_neighbor_run1=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    extend_neighbor_run1(i)=0;
    for j=1:len1
        extend_neighbor_run1(i)=run1(B(j))+extend_neighbor_run1(i);
    end
end

for i=1:N
    disp(i);
    j=1;%第一层
    nei0=[i];%第0层
    nei1=find(mixedsig(i,:));%算出节点i的直接邻居
    par=[];
    nei1=find(mixedsig(i,:));
    len1=length(nei1);
    par1=[];
    for j=1:len1
        par1(j)=core(i)*core(nei1(j));
    end
    
    par11=sum(par1);
    j=2;%第二层
    %     nei2=neighbor(mixedsig,nei1');
    % a=neighbor(mixedsig,nei1);
    a=[];
    for v=1:length(nei1)
        %     rr=mixedsig(nei1(v));
        k=find(mixedsig(nei1(v),:));
        a=union(k,a);
    end
    c=setdiff(a,nei0);%得到节点i的二跳邻居
    nei2=setdiff(c,nei1);%得到节点i的二跳邻居
    
    par2=[];
    for k=1:length(nei2)
        par2(k)=core(i)*core(nei2(k))/2^2;
        %         par2(k)=BBB(i,nei2(k))*(1/k1)^2+BBA(i,nei2(k))*(1/k1)^3;
    end
    par22=sum(par2);
    
    j=3;%第3层
    a=[];
    for v=1:length(nei2)
        %     rr=mixedsig(nei1(v));
        k=find(mixedsig(nei2(v),:));
        a=union(k,a);
    end
    c=setdiff(a,union(nei0,nei1));%得到节点i的3跳邻居
    nei3=setdiff(c,nei2);
    % BBA4=mixedsig^4;
    par3=[];
    for k=1:length(nei3)
        par3(k)=core(i)*core(nei3(k))/3^2;
    end
    par33=sum(par3);
    runk(i)=par11+par22+par33;
end
extend_neighbor_runk=[];
[r,l]=size(mixedsig);
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    extend_neighbor_runk(i)=0;
    for j=1:len1
        extend_neighbor_runk(i)=runk(B(j))+extend_neighbor_runk(i);
    end
end


% [Cm] = mixed_degree_decomposition(mixedsig, 0.3)
M_core=MDD(mixedsig);
core=core_numbers(sparse(mixedsig))';
avcore=sum(core)/N;
max_core=max(core);
h = h_index(mixedsig);
max_h=max(h);
av_h=sum(h)/N;
d = degreea( mixedsig );
max_d=max(d);
av_d=sum(d)/N;

% average_i=sir_main(mixedsig,2000,0.06) ;
% [B,index1]=sort(run1,'descend');
% cc1=index1(1,1:20);
% [B2,index2]=sort(rund,'descend');
% cc2=index2(1,1:20);
%  C11=setdiff(cc1,cc2)
% figure(1);
% ry=scatter(d,average_i,55,[101,147,74]/258,'s','filled');
% set(gca,'linewidth',1.8);
% box on;
% set(gca,'FontSize',15);
% xlabel('Degree (Enron)');
% ylabel('Φ(i)');
% text(0.1, 6.7, '(a)','FontSize',18);
% saveas(ry,'1Enron_degree.jpg');
% 
%  figure(2);
% ry=scatter(bet,average_i,55,[101,147,74]/258,'s','filled');
% set(gca,'linewidth',1.8);
% box on;
% set(gca,'FontSize',15);
% xlabel('Betweeness (Enron)');
% ylabel('Φ(i)');
% text(0.1, 6.7, '(b)','FontSize',18);
% saveas(ry,'2Enron_Betweeness.jpg');
% 
%  figure(3);
% ry=scatter(clo,average_i,55,[101,147,74]/258,'s','filled');
% set(gca,'linewidth',1.8);
% box on;
% set(gca,'FontSize',15);
% xlabel('Closeness (Enron)');
% ylabel('Φ(i)');
% text(0.1, 6.7, '(c)','FontSize',18);
% saveas(ry,'3Enron_Closeness.jpg');
% 
% 
%  figure(4);
% ry=scatter(rund,average_i,55,[101,147,74]/258,'s','filled');
% set(gca,'linewidth',1.8);
% box on;
% set(gca,'FontSize',15);
% xlabel('LGM (Enron)');
% ylabel('Φ(i)');
% text(0.1, 6.7, '(d)','FontSize',18);
% saveas(ry,'4Enron_LGM.jpg');
% 
% figure(5);
% ry=scatter(runk,average_i,55,[101,147,74]/258,'s','filled');
% set(gca,'linewidth',1.8);
% box on;
% set(gca,'FontSize',15);
% xlabel('G (Enron)');
% ylabel('Φ(i)');
% text(0.1, 6.7, '(e)','FontSize',18);
% saveas(ry,'5Enron_G.jpg');
% 
%  figure(6);
% ry=scatter(extend_neighbor_runk,average_i,55,[101,147,74]/258,'s','filled');
% set(gca,'linewidth',1.8);
% box on;
% set(gca,'FontSize',15);
% xlabel('G+ (Enron)');
% ylabel('Φ(i)');
% text(0.1, 6.7, '(f)','FontSize',18);
% saveas(ry,'6Enron_G+.jpg');
% 
%  figure(7);
% ry=scatter(runim,average_i,55,[101,147,74]/258,'s','filled');
% set(gca,'linewidth',1.8);
% box on;
% set(gca,'FontSize',15);
% xlabel('KSGC (Enron)');
% ylabel('Φ(i)');
% text(0.1, 6.7, '(g)','FontSize',18);
% saveas(ry,'7Enron_KSGC.jpg');
% 
%  figure(8);
% ry=scatter(run1,average_i,55,[101,147,74]/258,'s','filled');
% set(gca,'linewidth',1.8);
% box on;
% set(gca,'FontSize',15);
% xlabel('ISM (Enron)');
% ylabel('Φ(i)');
% text(0.1, 6.7, '(h)','FontSize',18);
% saveas(ry,'8Enron_ISM.jpg');
% 
%  figure(9);
% ry=scatter(extend_neighbor_run1,average_i,55,[101,147,74]/258,'s','filled');
% set(gca,'linewidth',1.8);
% box on;
% set(gca,'FontSize',15);
% xlabel('ISM+ (Enron)');
% ylabel('Φ(i)');
% text(0.1, 6.7, '(i)','FontSize',18);
% saveas(ry,'9Enron_ISM+.jpg');


%   xx=sir_main(mixedsig,500,0.02) ;
% %  xx=sir_mme(mixedsig,N,0.05,0,1000);
%  a00_clo=corr(xx',clo','type','Kendall');
% a00_bet=corr(xx',bet,'type','Kendall');
% a00_d=corr(xx',d','type','Kendall');
% a00_mdd=corr(xx',M_core','type','Kendall');
% a00_rund=corr(xx',rund','type','Kendall');
% a00_runim=corr(xx',runim','type','Kendall');
% a00_run1=corr(xx',run1','type','Kendall');
% a00_run1e=corr(xx',extend_neighbor_run1','type','Kendall');
% a00_runk=corr(xx',runk','type','Kendall');
% a00_runke=corr(xx',extend_neighbor_runk','type','Kendall');
% a00_clo=corr(average_i',clo','type','Kendall');
% a00_bet=corr(average_i',bet,'type','Kendall');
% a00_d=corr(average_i',d','type','Kendall');
% a00_mdd=corr(average_i',M_core','type','Kendall');
% a00_rund=corr(average_i',rund','type','Kendall');
% a00_runim=corr(average_i',runim','type','Kendall');
% a00_run1=corr(average_i',run1','type','Kendall');
% a00_run1e=corr(average_i',extend_neighbor_run1','type','Kendall');
% 
%  a00(i,1:12)=[a00_clo,a00_bet,a00_runk,a00_runke,a00_rund,a00_runim,a00_run1,a00_run1e,a00_run5,a00_run5e,a00_d,a00_neighbor_core];
 for i=4:15
   average_i{i}=sir_main(mixedsig,500,i*0.01) ;
a00_clo=corr(average_i{i}',clo','type','Kendall');
a00_bet=corr(average_i{i}',bet,'type','Kendall');
a00_d=corr(average_i{i}',d','type','Kendall');
a00_mdd=corr(average_i{i}',M_core','type','Kendall');
a00_runk=corr(average_i{i}',runk','type','Kendall');
a00_runke=corr(average_i{i}',extend_neighbor_runk','type','Kendall');
a00_rund=corr(average_i{i}',rund','type','Kendall');
% a00_runde=corr(average_i',extend_neighbor_rund','type','Kendall');
a00_runim=corr(average_i{i}',runim','type','Kendall');
% a00_runime=corr(average_i{i}',extend_neighbor_runim','type','Kendall');
a00_run1=corr(average_i{i}',run1','type','Kendall');
a00_run1e=corr(average_i{i}',extend_neighbor_run1','type','Kendall');
% a00(i,1:10)=[a00_clo,a00_bet,a00_d,a00_mdd,a00_rund,a00_runim,a00_run1,a00_run1e,a00_runk,a00_runke];
a00(i,1:8)=[a00_clo,a00_bet,a00_d,a00_mdd,a00_rund,a00_runim,a00_runk,a00_run1];
 end
 save('im_rt-twitter-copen.mat');
% a00=[];

% a00(:,11:14)=[];
% a00(:,12)=[];
% a00(:,8)=[];
% a00(:,6)=[];
% a00(:,3:4)=[];


