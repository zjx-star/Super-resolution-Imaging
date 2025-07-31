h=0.1;
i=1;
j=10;
R=2;
Out=[];
format long;
for hs=[h,h/2,h/2.^2,h/2.^3,h/2.^4]
beta=accu_beta_D(i,j,R,hs);
g=@(r)(bessely(i,beta*R).*besselj(i,beta.*r)-besselj(i,beta*R).*bessely(i,beta.*r)).^2.*(r);
C1=sqrt(2*pi*integral(g,R,R*(1+hs),'AbsTol',1e-8,'RelTol',1e-8));
C2=asymp_norm_psi_D(i,j,R,hs);
C3=norm_psi_D(i,j,R,hs);
Out=[Out;abs(C1-C2)/hs^(3/2),abs(C2-C3)/hs^(3/2),abs(C1-C3)/hs^(3/2)];
end
display(Out);