function d = degreea( mixedsig )
[l,n]=size(mixedsig)
for i=1:l
B=find(mixedsig(i,:));
d(i)=length(B);
end
end

