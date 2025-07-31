function val = H_mn_TE_3(x_3,r,theta,k,m,n,h)
%H_mn_TE Compute the value of H_mn_TE at certain positions
%Input:
%Precondition:
%Output:
%Postcondition:
m1=m;
m=abs(m);
s=s_mn_N(k,m,n,h);
beta=accu_beta_N(m,n,h);
lambda=eigenvalue_N(m,n,h);

g=@(r)((bessely(m-1,beta)-bessely(m+1,beta))./2.*besselj(m,beta.*r)-(besselj(m-1,beta)-besselj(m+1,beta))./2.*bessely(m,beta.*r)).^2;
C=sqrt(2*pi*integral(g,1,1+h,'AbsTol',1e-8,'RelTol',1e-8));
if(m==0)
    v=((bessely(m-1,beta)-bessely(m+1,beta))./2.*besselj(m,beta.*r)-(besselj(m-1,beta)-besselj(m+1,beta))./2.*bessely(m,beta.*r));
    val=1/k*-1i*2*cos(s*x_3)*lambda*v/C;
else
    v=((bessely(m-1,beta)-bessely(m+1,beta))./2.*besselj(m,beta.*r)-(besselj(m-1,beta)-besselj(m+1,beta))./2.*bessely(m,beta.*r))*exp(1i*m*theta);
    val=1/k*-1i*2*cos(s*x_3)*lambda*v/C;
end


end

