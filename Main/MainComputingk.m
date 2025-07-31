m=1;                                            %The index of the resonance
n=0;
N=5;                                            %The order of submatrices
R=0.6;
l=5;                                            %Thickness of the slab
h=0.1;                                          %The starting h
mode=0;                                         %0 imples even mode while 1 imples odd mode
diviser=2;
mu_ratio=1;
eps_ratio=25/9;
delta=1e-6;epsilon=1e-6;max1=10;

Out1=[];
i=1;
format long
for hs=[h/R h/R/diviser h/R/diviser^2]
   k0=kasymptotic(m,n,eps_ratio,mu_ratio,l,R,hs,mode)
   k1=real(k0);k2=k1+hs;k3=k1-hs; %Initial values
   %f2=@(x)(min(eig(LinearSystem(m,N,l,x,R,hs,mode))));
   f2=@(x)(min(eig(LinearSystem(m,N,l,x/sqrt(mu_ratio)/sqrt(eps_ratio),x,mu_ratio,R,hs,mode))));
   k=muller(f2,k1,k2,k3,delta,epsilon,max1);
   disp(k); 
   %Out
   Out1=[Out1;k,k0,k-k0,abs(k-k0)];
   i=i+1;
end
%%
%K=kasymptotic(m,[h h/diviser h/diviser^2 h/diviser^3 h/diviser^4 h/diviser^5]',l,VAlpha,MatAlpha,op);
%Klow=kasymptoticlow(1,[h h/diviser h/diviser^2 h/diviser^3]',1,0);
%format short
%disp([Out,K,Out-K,abs(Out-K)]);
disp([Out1]);


