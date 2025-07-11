
function kl=kliu(mixedsig,core)

D=mixedsig;
n=length(D);
% D(i,j)表示节点i和j的最短距离
for k=1:n
 for i=1:n
  for j=1:n
   if 0<D(i,k) & 0<D(k,j)
    if D(i,j)==0 & i~=j
      D(i,j)=D(i,k)+D(k,j);
    else 
      D(i,j)=min(D(i,j),D(i,k)+D(k,j));
    end
   end
  end
 end
end

max_core=max(core);
GG=find(core==max_core);
total_max_core=length(GG);
[r,l]=size(mixedsig);
kl=[];distance=[];

for i=1:r
    distance(i)=0;
for j=1:total_max_core
    distance(i)=D(i,GG(j))+distance(i)
end
end
for i=1:r
    kl(i)=(max_core-core(i)+1)*distance(i);
end
end