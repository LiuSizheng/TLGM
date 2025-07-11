%%
%%读取网络文件，找到最大连通分量
clear;
% 定义数据集文件名和对应的关键节点集合
datasets = {'1LFR2000-k5', 'LFR2000-10'};%'netscience','protain', 'facebook324', 'infectious', 'yeast1458', 'CA-GrQc', 'LFR2000-k15', 


initial_sets_Louvain = {
    %[4, 65, 67, 201, 21, 26, 40, 42, 170, 100, 214, 86, 113, 95, 128, 207, 270, 185, 334],%netscience
    %[1, 49, 763, 668, 107, 15, 127, 195, 1139, 179, 1346, 96, 261, 563, 583, 557, 266, 608, 139, 420, 151, 202, 705, 1331, 1082, 2207, 636, 591, 700, 784, 917, 1178, 1634],%protain
    %[107, 256, 182, 112, 232, 86, 284, 258, 15],%facebook
    %[51, 148, 122, 219, 292, 272],%infectious
    %[1312, 946, 638, 595, 27, 367, 1164, 819, 719, 72, 448, 567, 563, 345, 85, 1275, 223, 136, 333, 1130, 147, 129, 961, 458, 406, 397, 433, 1294, 454, 1315],%yeast
    %[51, 524, 1576, 3040, 891, 3470, 3823, 963, 221, 1100, 790, 3681, 2986, 2188, 607, 1949, 3013, 2396, 3134, 3275, 4126, 936, 291, 3380, 2669, 3244, 3262, 972, 3171, 2124, 3730, 2060, 1942, 823, 1889, 2926, 2033, 177, 3338, 3415, 4061],%CA-GrQc
    %[1977, 1909, 1824, 1956, 1979, 1992, 1941, 1965, 1853, 1980, 1988, 1964, 1779, 1968, 1976, 1647, 1934, 1948, 1952, 1707, 1975, 1815, 1843, 1756, 1999, 1715, 1900, 1932, 1990, 1991, 1919, 1680, 1930, 1971, 1884, 1962, 1902, 1966, 1933, 1727, 1989, 1914, 1921, 1997, 1978, 1951, 1982, 1957, 1622, 1949, 1985, 2000, 1938, 1927, 1863, 1825, 1724, 1833, 1954, 1891, 1963],%2000,15
    [1964, 1955, 1872, 1971, 1966, 1963, 1976, 1989, 1941, 1909, 1983, 1962, 1995, 1982, 1992, 1991, 1985, 1945, 1978, 1993, 1884, 1986, 1994, 2000, 1999, 1977, 1970, 1954, 1899, 1935, 1957, 1720, 1990, 1898, 1998, 1975, 1870, 1997, 1788, 1968, 1946, 1940, 1987, 1852, 1812, 1930, 1947, 1965, 1996, 1916],%2000,5
    [1849, 1846, 1945, 1999, 1947, 1869, 1920, 1925, 1974, 1977, 1995, 1903, 1986, 1978, 1940, 1872, 1953, 1989, 1902, 1993, 1901, 1909, 1900, 1790, 1856, 1988, 1976, 1950, 1960, 1908, 1980, 1998, 1984, 1814, 1997, 1857, 1973, 1914, 1975, 1756, 1948, 1982, 1659, 1881, 1961, 1842, 1904, 1853, 1919, 1892, 1877, 1983, 1815, 1911, 1992, 1956, 1985, 2000, 1962, 1991, 1757, 1759, 1935]%2000,10
};

initial_sets_LP = {
    %[4, 8, 10, 32, 26, 60, 40, 243, 275, 49, 70, 77, 88, 92, 95, 113, 118, 128, 143, 186, 211, 225, 236, 270, 335],%netscience
    %[608, 15, 17, 49, 26, 1510, 29, 2156, 70, 1085, 78, 110, 1971, 1295, 200, 263, 2207, 282, 304, 339, 367, 369, 763, 540, 426, 467, 446, 458, 459, 494, 1501, 638, 636, 933, 662, 729, 1182, 784, 861, 917, 993, 1016, 1061, 1152, 1202, 1245, 1449, 1843, 1582, 1634, 1644, 1658, 1681, 1700, 1817, 1854, 1998, 2113, 2200, 2282, 2258, 2344, 2477, 2489],%protain
    %[1, 256, 247, 104, 86, 258, 15, 232, 34, 284, 219, 274, 260, 308, 322],%facebook
    %[51, 61, 272],%infectious
    %[773, 2, 459, 595, 27, 367, 444, 98, 465, 719, 11, 329, 1105, 345, 114, 192, 20, 1212, 448, 24, 638, 36, 969, 567, 633, 542, 674, 147, 389, 264, 43, 45, 99, 458, 186, 52, 355, 1436, 267, 58, 906, 136, 63, 162, 72, 257, 70, 71, 451, 223, 683, 756, 436, 82, 433, 270, 85, 961, 715, 240, 819, 711, 96, 846, 406, 956, 134, 106, 131, 1099, 273, 1329, 129, 555, 271, 139, 141, 890, 530, 721, 255, 1113, 1235, 450, 755, 468, 229, 220, 920, 793, 397, 438, 361, 200, 335, 1020, 210, 212, 214, 647, 218, 844, 228, 318, 234, 235, 603, 466, 385, 1124, 1164, 319, 952, 287, 689, 306, 527, 309, 857, 550, 320, 708, 641, 364, 983, 946, 970, 1123, 383, 1025, 1192, 403, 412, 754, 420, 593, 1050, 778, 760, 492, 494, 736, 1125, 657, 512, 811, 1042, 559, 1074, 648, 579, 646, 650, 886, 680, 686, 697, 1315, 726, 1052, 1035, 790, 940, 1155, 1171, 1312, 1206, 1234, 1450, 1453],%yeast
    %[4017, 423, 1266, 4, 2420, 394, 2254, 3823, 3934, 117, 51, 2986, 2625, 1100, 2668, 790, 3897, 3681, 1290, 637, 2139, 26, 1949, 29, 3846, 3262, 2564, 674, 923, 2527, 39, 252, 3675, 42, 1910, 44, 2578, 891, 551, 1970, 524, 3380, 3013, 2669, 3549, 198, 1741, 1529, 3286, 1280, 426, 2217, 4158, 2253, 71, 2747, 75, 1778, 3767, 347, 139, 3205, 83, 2124, 1142, 1757, 88, 2643, 197, 3898, 146, 1210, 2854, 2700, 100, 2593, 1816, 516, 3218, 1207, 353, 1942, 1247, 3782, 2711, 2186, 3332, 127, 2061, 4073, 823, 2948, 2328, 745, 4085, 1316, 4126, 1617, 145, 451, 2409, 3481, 2579, 2861, 164, 385, 1997, 1576, 989, 4079, 2467, 177, 178, 350, 185, 1030, 1505, 3911, 3338, 291, 2122, 261, 4134, 1070, 2660, 1402, 3560, 943, 3530, 1768, 2569, 215, 216, 3666, 2993, 2990, 3485, 1839, 227, 607, 313, 4152, 803, 1844, 1930, 1625, 3237, 1999, 3171, 384, 269, 271, 3289, 4154, 812, 3002, 1639, 3577, 2331, 396, 3816, 1994, 733, 2772, 298, 1967, 808, 2647, 310, 1656, 2969, 322, 325, 326, 1628, 797, 1640, 3187, 657, 958, 2145, 3926, 3854, 1265, 1319, 1388, 3762, 363, 3517, 375, 376, 689, 3181, 433, 2033, 3590, 393, 2214, 399, 3077, 408, 2341, 415, 416, 2337, 596, 1192, 1140, 2316, 3818, 2237, 450, 439, 2365, 2940, 444, 453, 3064, 455, 3440, 1222, 1658, 735, 1774, 3155, 3658, 3869, 3371, 2725, 485, 3777, 1675, 492, 493, 3791, 691, 877, 2862, 569, 521, 3894, 2975, 963, 951, 2081, 3057, 1206, 2227, 2761, 2839, 3557, 1441, 1073, 1947, 2141, 579, 2654, 2714, 2536, 3500, 4061, 3964, 3929, 619, 1110, 1662, 3802, 2756, 925, 1704, 1592, 1197, 858, 3936, 682, 1152, 1338, 2987, 1058, 3699, 3441, 704, 3015, 1799, 3433, 722, 2461, 883, 3624, 2617, 841, 1961, 2026, 1055, 794, 800, 1318, 2416, 4054, 818, 821, 826, 4043, 2748, 832, 3089, 2434, 2099, 1755, 1166, 1889, 1583, 1973, 1706, 2239, 3985, 1870, 4008, 2633, 959, 3630, 956, 2361, 3613, 2089, 2544, 3083, 3085, 1010, 1663, 3080, 2635, 2924, 2638, 3880, 4089, 1328, 1765, 2411, 2172, 1833, 3374, 1155, 1188, 1194, 3111, 1219, 1225, 1226, 3076, 1269, 2868, 1293, 3122, 1306, 2749, 1957, 2475, 1771, 1908, 1376, 2389, 1654, 2458, 1745, 1477, 2769, 1809, 3769, 1650, 1563, 2575, 1584, 1602, 2865, 2082, 3615, 2225, 3809, 3178, 1940, 1752, 1880, 2218, 4002, 1827, 2750, 1859, 3303, 3362, 1926, 3030, 3665, 2717, 3098, 2498, 2667, 3377, 2179, 2639, 2343, 2811, 2425, 2503, 2583, 2732, 2792, 2933, 3065, 3143],%CA-GrQc
    %[1977, 1909, 1824, 1956, 1994, 1992, 1941, 1972, 1853, 1983, 1988, 1961, 1779, 1968, 1976, 1647, 1934, 1948, 1952, 1707, 1975, 1815, 1843, 1756, 1999, 1715, 1900, 1932, 1990, 1991, 1919, 1680, 1930, 1971, 1884, 1987, 1902, 1966, 1933, 1727, 1989, 1905, 1921, 1997, 1978, 1951, 1982, 1957, 1622, 1953, 1981, 1995, 1938, 1863, 1825, 1741, 1833, 1954, 1891, 1963],%2000,15
    [1964, 1762, 3, 1971, 1966, 1990, 1976, 1989, 1941, 1909, 1967, 1992, 1995, 1955, 1979, 1991, 1985, 1950, 1993, 1884, 1945, 1986, 1994, 1954, 33, 1970, 1999, 1977, 1874, 1899, 1257, 1957, 1898, 1978, 1962, 1700, 1691, 1998, 1618, 1981, 1870, 1914, 1997, 1825, 850, 1572, 1846, 1959, 85, 1872, 1668, 1358, 1940, 1987, 527, 106, 1852, 112, 1946, 534, 192, 575, 1812, 1930, 132, 1965, 143, 731, 857, 1972, 165, 1968, 394, 1836, 179, 189, 200, 1996, 1609, 210, 1916, 1129, 885, 1963, 223, 1761, 1935, 1944, 2000, 1947, 773, 270, 1709, 1583, 282, 768, 393, 610, 1600, 1586, 1699, 319, 320, 1261, 927, 1323, 1228, 1285, 447, 1249, 635, 1099, 505, 580, 1161, 642, 681, 930],%2000,5
    [1852, 1846, 1955, 1999, 1947, 1869, 1920, 1925, 1974, 1977, 1995, 1903, 1986, 1978, 1944, 1874, 1953, 1987, 1902, 1993, 1901, 1909, 1900, 1790, 1856, 1988, 1976, 1950, 1960, 1891, 1980, 1998, 1984, 1814, 1997, 1857, 1973, 1914, 1975, 1756, 1964, 1982, 1659, 1881, 1961, 1842, 1904, 1853, 1919, 1892, 1877, 1983, 1815, 1911, 1992, 1758, 1956, 1985, 1994, 1962, 1991, 1757, 1759, 1935]%2000,10
    
};


for dataset_idx = 1:length(datasets)
    % 清空上一次循环的变量
    clear mix degree pc core C ii totlenc num threshold initial_set kcore LI GAC GAC_Louvain GAC_LP
    filename = strcat(datasets{dataset_idx}, '.txt');
    net=load(filename); %读取网络文件
    adjMat=zeros(max(max(net))); %创建网络节点数量大小的邻接矩阵，初始为0
    %N=size(adjMat,1); %获取行数
    len=length(net); %获取net行数
    for i =1:len   %邻接赋值为1
        adjMat(net(i,1),net(i,2))=1;
        adjMat(net(i,2),net(i,1))=1;
    end
    %spMat=sparse(adjMat);%稀疏矩阵
    %[a b]=components(spMat);%找到连通分量 a-连通分量序号 b-分量节点数
    [B]=largestcomponent(adjMat); %B返回最大连通分量包含的节点序号
    adjMat=adjMat(B,B); %得到最大连通分量的邻接矩阵
    mix=adjMat;
    N=size(mix,1);  %节点数
    degree=sum(mix,2)'; %按行求和，得节点的度
    avgdeg=sum(degree)/N;  %平均度
    pc=zeros(1,N);
    for i=1:N
        B1=find(mix(i,:));%和第i个节点相连的节点索引
        len1=length(B1);%len1表示i的邻居数
        ddin(i)=0;
        for j=1:len1
            B2=find(mix(B1(j),:));%i的邻居之一j的邻居
            jiaoji=intersect(B1,B2);%i与j的邻居交集
            ddin(i)=length(jiaoji)+ddin(i);
        end
        pc(i)=degree(i)*sum(ddin(i));
    end
    core=pc;
    C={};
    ii=1;
    totlenc=0;
    num=0;
    threshold=0;
    while threshold<=0.6  %达到网络节点总数的60%停止
        %[a,b]=sort(core,'descend'); %a-降序序列 b-节点索引
        [pcmax,index]=max(core);
        biggest=find(core==pcmax);%找到pcmax的节点
        aaa=biggest(find(degree(biggest)==max(degree(biggest))));%找到pc值最大的里面度最大的
        len0=length(aaa);
        initial_set(ii)=aaa(ceil(rand(1)*len0));%随机取得最大pc值中的一个作为起始节点%%%%
        num=num+1;
        if num/N<=0.6
            %initial_set(ii)=aaa(1);
            C{ii}=initial_set(ii);
            neighbor_of_i=find(mix(initial_set(ii),:))%初始节点的邻居都找出来
            zz=find(pc(neighbor_of_i)>=pc(initial_set(ii))/2);%找到pc值大于1/2pcmax的节点
            pc_of_zz=pc(neighbor_of_i(zz));
            togther0=[pc_of_zz;neighbor_of_i(zz)]';
            togtherr0=sortrows(togther0)';
            for j=length(zz):-1:1
                C{ii}=[C{ii},togtherr0(2,j)];
                num=num+1;
                if num/N>0.6
                    C{ii} = C{ii}(1:end-1);
                    break
                end
            end
            % C{ii}=neighbor_of_i(zz);%找到节点i的内部链接的各个节点,组成子图C
            % num=num+length(C{ii});
            %C{ii}=[C{ii},initial_set(ii)];
            if num/N<=0.6
                first_neighbor_of_C=neighbor_of_node_set(mix,C{ii});
                degree_of_first_neighbor_of_C=degree(first_neighbor_of_C);
                togther=[degree_of_first_neighbor_of_C;first_neighbor_of_C]';
                togtherr=sortrows(togther)';
    
                %for j=length(first_neighbor_of_C):-1:1
                for j=1:length(first_neighbor_of_C)
                    ddd=neighbor(mix,togtherr(2,j));
                    jiaoji=intersect(ddd,C{ii});
                    chaji=setdiff(ddd,jiaoji);
                    if length(jiaoji)>length(chaji)
                        C{ii}=[C{ii},togtherr(2,j)];
                        num=num+1;
                        if num/N>0.6
                            C{ii} = C{ii}(1:end-1);
                            break
                        end
                    end
                end
            else
                break
            end
            if num/N<=0.6
                second_neighbor_of_C=neighbor_of_node_set(mix,C{ii});
                degree_of_second_neighbor_of_C=degree(second_neighbor_of_C);
                togther=[degree_of_second_neighbor_of_C;second_neighbor_of_C]';
                togtherr=sortrows(togther)';
                %for j=length(second_neighbor_of_C):-1:1
                for j=1:length(second_neighbor_of_C)
                    ddd=neighbor(mix,togtherr(2,j));
                    jiaoji=intersect(ddd,C{ii});
                    chaji=setdiff(ddd,jiaoji);
                    if length(jiaoji)>length(chaji)
                        C{ii}=[C{ii},togtherr(2,j)];
                        num=num+1;
                        if num/N>0.6
                            C{ii} = C{ii}(1:end-1);
                            break
                        end
                    end
                end
            else
                break
            end
            if num/N<=0.6
                third_neighbor_of_C=neighbor_of_node_set(mix,C{ii});
                degree_of_third_neighbor_of_C=degree(third_neighbor_of_C);
                togther=[degree_of_third_neighbor_of_C;third_neighbor_of_C]';
                togtherr=sortrows(togther)';
                %for j=length(third_neighbor_of_C):-1:1
                for j=1:length(third_neighbor_of_C)
                    ddd=neighbor(mix,togtherr(2,j));
                    jiaoji=intersect(ddd,C{ii});
                    chaji=setdiff(ddd,jiaoji);
                    if length(jiaoji)>length(chaji)
                        C{ii}=[C{ii},togtherr(2,j)];
                        num=num+1;
                        if num/N>0.6
                            C{ii} = C{ii}(1:end-1);
                            break
                        end
                    end
                end
            else
                break
            end
            %C_number_of_initial_i(i)=length(C{i});
            linshi3=C{ii};
            lenc=length(linshi3);
            totlenc=totlenc+lenc;
            linshi4=linshi3(find(degree(linshi3)==max(degree(linshi3))));
            initial_set(ii)=linshi4(ceil(rand(1)*length(linshi4)));
            %initial_set(ii)=linshi4(1);
            mix(C{ii},:)=0;
            mix(:,C{ii})=0;
            degree=sum(mix);
            MN=sparse(mix);
            for i=1:N
                B1=find(mix(i,:));%和第i个节点相连的节点索引
                len1=length(B1);%len1表示i的邻居数
                ddin(i)=0;
                for j=1:len1
                    B2=find(mix(B1(j),:));%i的邻居之一j的邻居
                    jiaoji=intersect(B1,B2);%i与j的邻居交集
                    ddin(i)=length(jiaoji)+ddin(i);
                end
                core(i)=degree(i)*sum(ddin(i));
            end
            ii=ii+1;
            threshold=totlenc/N;    
        else
            initial_set=initial_set(1:end-1);
            break
        end
    end
    %%节点局部影响计算
    kcore=core_numbers(sparse(adjMat))';
    LI=[];
    for i=1:N
        B3=find(adjMat(i,:));%和第i个节点相连的节点索引
        len2=length(B3);%len3表示i的邻居数
        ks=0;
        for j=1:len2
            ks=ks+kcore(B3(j));
        end
        LI(end+1)=ks;
    end

    shortest_distances = dijkstra_for_targets(adjMat,initial_set);%找到全局关键节点到其他节点的最短距离
    sd=shortest_distances';
    for i=1:N
        GAC(i)=0;
        for j=1:length(initial_set)
            GAC(i)=LI(i)*LI(initial_set(j))/(2^sd(i,j))+GAC(i);
        end
    end

    initial_set_Louvain = initial_sets_Louvain{dataset_idx};
    shortest_distances = dijkstra_for_targets(adjMat,initial_set_Louvain);%找到全局关键节点到其他节点的最短距离
    sd=shortest_distances';
    GAC_Louvain = zeros(1,N);
    for i=1:N
        for j=1:length(initial_set_Louvain)
            GAC_Louvain(i)=LI(i)*LI(initial_set_Louvain(j))/(2^sd(i,j))+GAC_Louvain(i);
        end
    end

    initial_set_LP=initial_sets_LP{dataset_idx};
    shortest_distances = dijkstra_for_targets(adjMat,initial_set_LP);%找到全局关键节点到其他节点的最短距离
    sd=shortest_distances';
    GAC_LP = zeros(1,N);
    for i=1:N
        for j=1:length(initial_set_LP)
            GAC_LP(i)=LI(i)*LI(initial_set_LP(j))/(2^sd(i,j))+GAC_LP(i);
        end
    end
    % 保存结果到 .mat 文件
    filename = strcat(datasets{dataset_idx}, '_GAC_Louvain_LP.mat');
    savGACe(filename, 'initial_set','initial_set_Louvain','initial_set_LP','GAC','GAC_Louvain', 'GAC_LP');
end





