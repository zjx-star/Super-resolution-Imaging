h=0.2;
m=0;
n=1;
R=2;
Out1=[];
Out2=[];
format long;
for hs=[h,h/2,h/2.^2,h/2.^3,h/2.^4]
beta=accu_beta_N(m,n,R,hs);
g=@(r)((bessely(m-1,R*beta)-bessely(m+1,R*beta))/2*besselj(m,beta.*r)-(besselj(m-1,R*beta)-besselj(m+1,R*beta))/2*bessely(m,beta.*r)).^2;
C1=sqrt(2*pi*R*integral(g,R,R*(1+hs),'AbsTol',1e-8,'RelTol',1e-8));
C2=asymp_norm_psi_N(m,n,R,hs);
C3=norm_psi_N(m,n,R,hs);
Out1=[Out1;abs(C1-C2)/hs^(3/2),abs(C2-C3)/hs^(3/2)];
end
display(Out1);

m=1;
n=0;
for hs=[h,h/2,h/2.^2,h/2.^3,h/2.^4]
beta=accu_beta_N(m,n,R,hs);
g=@(r)((bessely(m-1,R*beta)-bessely(m+1,R*beta))/2*besselj(m,beta.*r)-(besselj(m-1,R*beta)-besselj(m+1,R*beta))/2*bessely(m,beta.*r)).^2;
C1=sqrt(2*pi*R*integral(g,R,R*(1+hs),'AbsTol',1e-8,'RelTol',1e-8));
C2=asymp_norm_psi_N(m,n,R,hs);
C3=norm_psi_N(m,n,R,hs);
Out2=[Out2;abs(C1-C2)/hs^(3/2),abs(C2-C3)/hs^(3/2)];
end
display(Out2);