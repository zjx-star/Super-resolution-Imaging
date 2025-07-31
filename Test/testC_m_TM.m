h=0.1;
m=-1;
N=5;
R1=2;
R2=2;
mode=0;
k=sqrt(1);
k1=sqrt(2);
Out=[];
format long;
for hs=[h,h/4,h/4.^2]
C_m_TM1=asymp_C_TM_R(m,N,k,k1,R1,R2,hs,mode);
C_m_TM2=C_m_TM(m,N,1,k,k1,R1,R2,hs,hs,mode);
Out=[Out;((C_m_TM1-C_m_TM2).*((C_m_TM1-C_m_TM2).')').^(1/2)];
end
display(Out);