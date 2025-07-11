function hefei1= GAC( mix,N,r,ddin,threshold,biggest_core,DegreeAll,totlenc,mixedsig,neighbor_Core )
i=1;
while threshold<=r
%for i=1:K
    
    aaa=biggest_core(find(DegreeAll(biggest_core)==max(DegreeAll(biggest_core))));
    len0=length(aaa);
    initial_set(i)=aaa(ceil(rand(1)*len0));%随机取得最大核数中的一个作为起始节点%%%%
    neighbor_of_i=find(mix(initial_set(i),:))%初始节点的邻居都找出来
    zz=find(ddin(neighbor_of_i)>=ddin(initial_set(i))/2);
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
    lenc=length(linshi3);
    totlenc=totlenc+lenc;
    linshi4=linshi3(find(DegreeAll(linshi3)==max(DegreeAll(linshi3))));
    initial_set(i)=linshi4(ceil(rand(1)*length(linshi4)));
    
    mix(C{i},:)=0;
    mix(:,C{i})=0;
    [a,b]=sort(ddin,'descend');
    MN=sparse(mix);
    core=iniddin(mix,N);
    %core=core_numbers(MN)';
    [a,b]=sort(core,'descend');
    [max_core,indexm]=max(core);%找出第i轮中得分最高的那个值
    biggest_core=find(core==max_core);
    len0=length(biggest_core);
    DegreeAll=sum(mix);
    i=i+1;  
    threshold=totlenc/N;
end
lenghefei=length(initial_set);
core=core_numbers(sparse(mixedsig))';
 h = (h_index(mixedsig))';
 neighbor_h=[];
for i=1:r
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    neighbor_h(i)=0;
    for j=1:len1
        neighbor_h(i)=core(B(j))+neighbor_h(i);
    end
end
shortest_distances = dijkstra_for_targets(mixedsig, initial_set);
sd=shortest_distances';
for i=1:N
    hefei1(i)=0;
    for j=1:lenghefei
    %DD = dijkstra(mixedsig,initial_set(j));
    %sssss=core(i);
    %dfs3e=core(initial_set(j));
    %    hefei1(i)=neighbor_Core(i)*neighbor_Core(initial_set(j))/(2^result(i,initial_set(j)))+hefei1(i);
            hefei1(i)=neighbor_Core(i)*neighbor_Core(initial_set(j))/(2^sd(i,j))+hefei1(i);
                    %    hefei1(i)=(neighbor_h(i)+neighbor_Core(i))/1*(neighbor_h(initial_set(j))*neighbor_Core(initial_set(j))/1/(2^sd(i,j))+hefei1(i));

               % hefei1(i)=exp(neighbor_Core(i))*exp(neighbor_Core(initial_set(j)))/(2^sd(i,j))+hefei1(i);

    end
end
end