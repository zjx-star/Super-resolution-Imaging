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
for hs=[h,h/4,h/4.^2]
    C_m_TE1=asymp_C_TE_R(m,N,k,k1,R1,R2,hs,mode);
    C_m_TE2=C_m_TE(m,N,1,k,k1,R1,R2,hs,hs,mode);
    Out=[Out;((C_m_TE1-C_m_TE2).*((C_m_TE1-C_m_TE2).')').^(1/2)];
end
display(Out);