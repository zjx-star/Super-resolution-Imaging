h=0.2;
m=1;
mode=0;
R1=2;
R2=2;
k=sqrt(1);
k1=sqrt(2);
Out=[];
format long;
for hs=[h,h/2,h/2.^2,h/2.^3]
A_mm1=asymp_A_mm_R(m,k,k1,R1,R2,hs,mode);
A_mm2=A_mm(m,1,k,k1,R1,R2,hs,hs,mode);
%A_mm2=asymp_A_mm(m,k,hs,mode)
delta=A_mm1-A_mm2;
Out=[Out;((delta(1:end)).*((delta(1:end)).')').^(1/2)];
end
display(Out);