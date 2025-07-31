h=0.2;
m=0;
n=1;
R=2;
Out1=[];
Out2=[];
Out3=[];
Out4=[];
format long;
psi_N=@(r,theta,h)exp(m*1i*theta)./sqrt(pi.*R.*r.*h).*cos(n.*pi./R.*(r-R)./h);
d_psi_N=@(r,theta,h)-n*pi^(1/2)/h^(3/2).*exp(m*1i*theta)./R.^2.*sin(n.*pi./R.*(r-R)./h);
for hs=[h,h/4,h/4.^2,h/4.^3]
beta=accu_beta_N(m,n,R,hs);
g=@(r)((bessely(m-1,R*beta)-bessely(m+1,R*beta))/2*besselj(m,beta.*r)-(besselj(m-1,R*beta)-besselj(m+1,R*beta))/2*bessely(m,beta.*r)).^2;
C1=sqrt(2*pi*R*integral(g,R,R*(1+hs),'AbsTol',1e-8,'RelTol',1e-8));
psi_N2=@(r,theta,h)((bessely(m-1,R*beta)-bessely(m+1,R*beta))/2*besselj(m,beta.*r)-(besselj(m-1,R*beta)-besselj(m+1,R*beta))/2*bessely(m,beta.*r))./C1*exp(1i*m*theta);
d_psi_N2=@(r,theta,h)beta.*((bessely(m-1,R*beta)-bessely(m+1,R*beta))/2*(besselj(m-1,beta.*r)-besselj(m+1,beta.*r))/2-(besselj(m-1,R*beta)-besselj(m+1,R*beta))/2*(bessely(m-1,beta.*r)-bessely(m+1,beta.*r))/2)./C1*exp(1i*m*theta);
Out1=[Out1,norm(psi_N(R*(1:0.1*hs:1+hs),0,hs)-psi_N2(R*(1:0.1*hs:1+hs),0,hs),2)];
Out2=[Out2,norm((d_psi_N(R*(1:0.1*hs:1+hs),0,hs)-d_psi_N2(R*(1:0.1*hs:1+hs),0,hs)),2)];
end
display(Out1);
display(Out2);

m=1;
n=0;
psi_N=@(r,theta,h)exp(m*1i*theta)./sqrt(pi.*h.*(h+2))./R;
d_psi_N=@(r,theta,h)0;
for hs=[h,h/4,h/4.^2,h/4.^3]
beta=accu_beta_N(m,n,R,hs);
g=@(r)((bessely(m-1,R*beta)-bessely(m+1,R*beta))/2*besselj(m,beta.*r)-(besselj(m-1,R*beta)-besselj(m+1,R*beta))/2*bessely(m,beta.*r)).^2;
C1=sqrt(2*pi*R*integral(g,R,R*(1+hs),'AbsTol',1e-8,'RelTol',1e-8));
psi_N2=@(r,theta,h)((bessely(m-1,R*beta)-bessely(m+1,R*beta))/2*besselj(m,beta.*r)-(besselj(m-1,R*beta)-besselj(m+1,R*beta))/2*bessely(m,beta.*r))./C1*exp(1i*m*theta);
d_psi_N2=@(r,theta,h)beta.*((bessely(m-1,R*beta)-bessely(m+1,R*beta))/2*(besselj(m-1,beta.*r)-besselj(m+1,beta.*r))/2-(besselj(m-1,R*beta)-besselj(m+1,R*beta))/2*(bessely(m-1,beta.*r)-bessely(m+1,beta.*r))/2)./C1*exp(1i*m*theta);
Out3=[Out3,norm(psi_N(R*(1:0.1*hs:1+hs),0,hs)-psi_N2(R*(1:0.1*hs:1+hs),0,hs),2)];
Out4=[Out4,norm((d_psi_N(R*(1:0.1*hs:1+hs),0,hs)-d_psi_N2(R*(1:0.1*hs:1+hs),0,hs)),2)];
end
display(Out3);
display(Out4);

