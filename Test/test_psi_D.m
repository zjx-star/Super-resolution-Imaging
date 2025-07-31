h=0.2;
i=0;
j=1;
R=2;
Out1=[];
Out2=[];
format long;
psi_D=@(r,theta,h)exp(i*1i*theta)./sqrt(pi.*R.*r.*h).*sin(j.*pi./R.*(R-r)./h);
d_psi_D=@(r,theta,h)-j*pi^(1/2)/h^(3/2).*exp(i*1i*theta)./R.^2.*cos(j.*pi./R.*(R-r)./h);
for hs=[h,h/4,h/4.^2,h/4.^3]
beta=accu_beta_D(i,j,R,hs);
g=@(r)(bessely(i,beta*R).*besselj(i,beta.*r)-besselj(i,beta*R).*bessely(i,beta.*r)).^2;
C1=sqrt(2*pi*R*integral(g,R,R*(1+hs),'AbsTol',1e-8,'RelTol',1e-8));
psi_D2=@(r,theta,h)(bessely(i,beta*R).*besselj(i,beta.*r)-besselj(i,beta*R).*bessely(i,beta.*r))./C1*exp(1i*i*theta);
d_psi_D2=@(r,theta,h)beta.*(bessely(i,beta*R).*(besselj(i-1,beta.*(r))-besselj(i+1,beta.*(r)))./2-besselj(i,beta*R).*(bessely(i-1,beta.*r)-bessely(i+1,beta.*(r)))./2)./C1*exp(1i*i*theta);
Out1=[Out1,norm(psi_D(R*(1:0.1*hs:1+hs),0,hs)-psi_D2(R*(1:0.1*hs:1+hs),0,hs),2)];
Out2=[Out2,norm((d_psi_D(R*(1:0.1*hs:1+hs),0,hs)-d_psi_D2(R*(1:0.1*hs:1+hs),0,hs)),2)];
end
display(Out1);
display(Out2);