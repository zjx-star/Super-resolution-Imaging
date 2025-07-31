h=0.2;
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
B_m_TM_TM1=asymp_B_TM_TM(m,N,k,k1,1,R1,R2,hs,mode);
B_m_TM_TM2=B_m_TM_TM(m,N,1,k,k1,R1,R2,hs,hs,mode);
delta=B_m_TM_TM1-B_m_TM_TM2;
Out=[Out;((delta(1:end)).*((delta(1:end)).')').^(1/2)];
end
display(Out);