m=1;
n=0;
h=0.1;
x3=0.5;
theta=pi/3;
Out=[];

E1=@(h)2i*m*exp(1i*m*theta)*cos(theta)/sqrt(2*pi*h);
H3=@(h)-2*m^2*1i*exp(1i*m*theta)/sqrt(2*pi*h);

for hs=[h,h/4,h/4.^2,h/4.^3]
k=m;
r=1+hs/2;
E1_asymp=E1(hs);
E_accu=E_mn_TE(x3,r,theta,k,m,n,hs);
E1_accu=E_accu(1);
delta1=E1_asymp-E1_accu;
H3_asymp=H3(hs);
H_accu=H_mn_TE_3(x3,r,theta,k,m,n,hs);
H3_accu=H_accu(1);
delta2=H3_asymp-H3_accu;
ratio1=E1_accu/H3_accu;
ratio2=E_accu(2)/H3_accu;
Out=[Out;((delta1(1:end)).*((delta1(1:end)).')').^(1/2),((delta2(1:end)).*((delta2(1:end)).')').^(1/2),ratio1,ratio2];
end
display(Out);


