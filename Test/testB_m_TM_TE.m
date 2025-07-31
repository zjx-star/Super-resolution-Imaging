h=0.1;
m=1;
N=3;
R1=1;
R2=1;
mode=0;
k=sqrt(1);
k1=sqrt(2);
Out=[];
format long;
for hs=[h,h/2,h/2.^2]
B_m_TM_TE1=zeros(N,N);
B_m_TM_TE2=B_m_TM_TE(m,N,1,k,k1,R1,R2,hs,hs,mode);
delta=B_m_TM_TE1-B_m_TM_TE2;
Out=[Out;((delta(1:end)).*((delta(1:end)).')').^(1/2)];
end
display(Out);