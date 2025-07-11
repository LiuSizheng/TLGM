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
%          runim(i)=par11+par22+par33;
   KSGC(i)=par11;
end
h = h_index(mixedsig);
c=pagerank(mixedsig,N);
for i=1:15
%average_i=sir_main(mixedsig,1000,0.05) ;
%average_i=average_iall(1) ;
a00_degree(i)=corr(average_iall{i}',d','type','Kendall');
a00_h(i)=corr(average_iall{i}',h,'type','Kendall');
a00_pagerank(i)=corr(average_iall{i}',c,'type','Kendall');
a00_M_core(i)=corr(average_iall{i}',M_core','type','Kendall');
a00_bet(i)=corr(average_iall{i}',bet,'type','Kendall');
a00_clo(i)=corr(average_iall{i}',clo','type','Kendall');
%a00_densityim=corr(average_i',runkim','type','Kendall');
a00_hefei1(i)=corr(average_iall{i}',hefei1','type','Kendall');
a00_extend_hefei(i)=corr(average_iall{i}',extend_hefei','type','Kendall');
a0_core(i)=corr(average_iall{i}',core','type','Kendall');
a00_extend_neighbor_core(i)=corr(average_iall{i}',extend_neighbor_core','type','Kendall')
a00_neighbor_core(i)=corr(average_iall{i}',neighbor_Core','type','Kendall');
a00_runk(i)=corr(average_iall{i}',runk','type','Kendall');
a00_KSGC(i)=corr(average_iall{i}',KSGC','type','Kendall');
end
all_average_sir=[a00_hefei1;a00_M_core;a00_bet;a00_clo;a00_extend_hefei;a00_neighbor_core;a00_degree;a0_core;a00_extend_neighbor_core;a00_runk;a00_KSGC];
all_average_sir2=[a00_hefei1;a00_h;a00_pagerank;a00_M_core;a00_neighbor_core;a00_degree;a0_core;a00_KSGC];
all_average_sir3=all_average_sir2';