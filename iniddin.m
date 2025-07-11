function ddin=iniddin(mixedsig,N)
for i=1:N
    B=find(mixedsig(i,:));%find函数用来返回矩阵mixedsig（i,:）中非零元素所在位置
    len1=length(B);%len1表示数组B的长度
    ddin(i)=0;
    for j=1:len1
         B2=find(mixedsig(B(j),:));
          jiaoji=intersect(B,B2);
        ddin(i)=degree(i,mixedsig)*length(jiaoji)+ddin(i);
    end
end
end