h=0.2;
m=1;
N=3;
R1=1;
R2=1;
mode=1;
k=sqrt(1);
k1=sqrt(2);
Out=[];
format long;
for hs=[h,h/4,h/4.^2]
R_m_TM1=asymp_R_TM_R(m,N,k,k1,R1,R2,hs,mode);
R_m_TM2=R_m_TM(m,N,1,k,k1,R1,R2,hs,hs,mode);
Out=[Out;((R_m_TM1-R_m_TM2).*((R_m_TM1-R_m_TM2).')').^(1/2)];
end
display(Out);