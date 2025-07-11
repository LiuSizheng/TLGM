function nc=NC(mixedsig)
mix=mixedsig;
N = size(mix, 1); % 当前图的节点数
d = degreea(mix);
nc=zeros(N,1);
for i=1:N
    neibor=find(mix(i,:));%和第i个节点相连的节点索引
    len2=length(neibor);%len3表示i的邻居数
    for j=1:len2
        nc(i)=nc(i)+d(neibor(j));
    end
end