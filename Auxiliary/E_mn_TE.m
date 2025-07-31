function mat = E_mn_TE(x_3,r,theta,k,m,n,h)
%E_mn_TE Compute the value of E_mn_TE at certain positions
%Input:
%Precondition:
%Output:
%Postcondition:
m1=m;
m=abs(m);
s=s_mn_N(k,m,n,h);
beta=accu_beta_N(m,n,h);

g=@(r)((bessely(m-1,beta)-bessely(m+1,beta))./2.*besselj(m,beta.*r)-(besselj(m-1,beta)-besselj(m+1,beta))./2.*bessely(m,beta.*r)).^2;
C=sqrt(2*pi*integral(g,1,1+h,'AbsTol',1e-8,'RelTol',1e-8));
if(m==0)
    v=beta*((bessely(m-1,beta)-bessely(m+1,beta))./2.*(besselj(m-1,beta.*r)-besselj(m+1,beta.*r))./2-(besselj(m-1,beta)-besselj(m+1,beta))./2.*(bessely(m-1,beta.*r)-bessely(m+1,beta.*r))./2);
    mat=2*cos(s*x_3)*v/C*[sin(theta),-cos(theta),0];
else
    v1=beta*((bessely(m-1,beta)-bessely(m+1,beta))./2.*(besselj(m-1,beta.*r)-besselj(m+1,beta.*r))./2-(besselj(m-1,beta)-besselj(m+1,beta))./2.*(bessely(m-1,beta.*r)-bessely(m+1,beta.*r))./2)*exp(1i*m*theta);
    v2=(1i*m)*((bessely(m-1,beta)-bessely(m+1,beta))./2.*besselj(m,beta.*r)-(besselj(m-1,beta)-besselj(m+1,beta))./2.*bessely(m,beta.*r))*exp(1i*m*theta);
    mat=2*cos(s*x_3)*(v1/C*[sin(theta),-cos(theta),0]+v2/C*[cos(theta)/r,sin(theta)/r,0]);
end


end

