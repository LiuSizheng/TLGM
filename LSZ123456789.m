clear;
% 加载之前保存的数组

load('LFR2000-k15_average_sir1_15.mat');
load('LFR2000-k15.mat');
for i=1:15
    a00_degree(i)=corr(average_iall{i}',d','type','Kendall');%度中心性
    a00_h_index(i)=corr(average_iall{i}',h,'type','Kendall');%H指数
    a00_pagerank(i)=corr(average_iall{i}',p,'type','Kendall');%PageRank中心性
    a00_MDD_core(i)=corr(average_iall{i}',MDD_core','type','Kendall');%混合度分解法MDD
    a00_GAC(i)=corr(average_iall{i}',GAC','type','Kendall');
    a00_kcore(i)=corr(average_iall{i}',kcore','type','Kendall');%k壳分解法
    a00_neighbor_core(i)=corr(average_iall{i}',neighbor_Core','type','Kendall');%邻域核心度中心性（NCC）
    a00_KSGC(i)=corr(average_iall{i}',KSGC','type','Kendall');%基于 K 核分解法的引力模型改进方法KSGC
end

save('LFR2000-k15_sir1_15_Kendall.mat', 'a00_degree', 'a00_h_index','a00_pagerank','a00_MDD_core','a00_GAC','a00_kcore','a00_neighbor_core','a00_KSGC');
