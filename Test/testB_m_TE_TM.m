h=0.1;
m=1;
N=3;
mode=0;
k=sqrt(1);
k1=sqrt(2);
Out=[];
format long;
for hs=[h,h/2,h/2.^2]
B_m_TE_TM1=zeros(N,N);
B_m_TE_TM2=B_m_TE_TM(m,N,1,k,k1,hs,mode);
delta=B_m_TE_TM1-B_m_TE_TM2;
Out=[Out;((delta(1:end)).*((delta(1:end)).')').^(1/2)];
end
display(Out);