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
for hs=[h,h/4,h/4.^2,h/4.^3]
R_m_TE1=asymp_R_TE_R(m,N,k,k1,R1,R2,hs,mode);
R_m_TE2=R_m_TE(m,N,1,k,k1,R1,R2,hs,hs,mode);
delta=R_m_TE1-R_m_TE2;
Out=[Out;((delta).*((delta).')').^(1/2)];
end
display(Out);