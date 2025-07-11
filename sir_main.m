function average_i=sir_main(mixedsig,n,irate) 

  A=mixedsig;
  N=size(A,1); 
  b=sum(A,2);
  k1=sum(b)/N;
  k2=sum(b'*b)/N; 
  %temp=k1/k2;
%irate=k1/k2;
%irate=0.15;
rrate=1;
for k=1:n
    disp(k);
for r=1:N
[r_result,i_result]=BFSspreading(N,A,r,irate,rrate);
infection_result(k,r)=i_result;
recover_result(k,r)=r_result;
end
end

NN=sum(recover_result);
average_i=NN/n;
 [rr,l]=sort(average_i) ;
end
    

