function vec = TE_RHS(m,N,l,k,R1,h1,mode,d1)
%TE_RHS This function generates the right hand side of the linear system 
%describing the excitation of TE mode when a normal incident plane wave is 
%provided.
%
%Input:
%Precondition:
%Output:
%Postcondition:

m1=m;
m=abs(m1);
vec1=zeros((2*N+1),1);
if(d1==0)
    parfor n=0:N
        beta1=accu_beta_N(m,n,R1,h1);
        lambda1=eigenvalue_N(m,n,R1,h1);
        s1=s_mn_N(k,m,n,R1,h1);
        f1=@(s1,theta)R1.*(1+s1.*h1).*...
            ((bessely(m-1,R1.*beta1)-bessely(m+1,R1.*beta1))./2.*besselj(m,beta1.*R1.*(1+s1.*h1))-(besselj(m-1,R1.*beta1)-besselj(m+1,R1.*beta1))./2.*bessely(m,beta1.*R1.*(1+s1.*h1))).*...
            R1.*(1+s1.*h1).*cos(-m1.*theta).*cos(theta);
        g=@(r)((bessely(m-1,R1.*beta1)-bessely(m+1,R1.*beta1))./2.*besselj(m,beta1.*r)-(besselj(m-1,R1.*beta1)-besselj(m+1,R1.*beta1))./2.*bessely(m,beta1.*r)).^2;
        C=sqrt(2*pi*integral(g,R1,R1*(1+h1),'AbsTol',1e-8,'RelTol',1e-8));
        vec1(n+1)=-k/(2i)/s1*lambda1.^(3/4)*exp(-1i*k*l/2)*R1*h1*integral2(f1,0,1,0,2*pi)/C;              
    end
    vec1(1)=vec1(1)*s_mn_N(k,m,0,R1,h1);
    if(mode==1)
        vec1=1i*vec1;
    end
    vec=[vec1];
else
    d3=sqrt(1-d1.^2);
    for n=0:N
        beta1=accu_beta_N(m,n,R1,h1);
        lambda1=eigenvalue_N(m,n,R1,h1);
        s1=s_mn_N(k,m,n,R1,h1);
        f1=@(s1,theta)R1.*(1+s1.*h1).*...
            ((bessely(m-1,R1.*beta1)-bessely(m+1,R1.*beta1))./2.*besselj(m,beta1.*R1.*(1+s1.*h1))-(besselj(m-1,R1.*beta1)-besselj(m+1,R1.*beta1))./2.*bessely(m,beta1.*R1.*(1+s1.*h1))).*...
            (exp(1i.*k.*d1.*R1.*(1+s1.*h1).*cos(theta))).*cos(-m1.*theta);
        g=@(r)((bessely(m-1,R1.*beta1)-bessely(m+1,R1.*beta1))./2.*besselj(m,beta1.*r)-(besselj(m-1,R1.*beta1)-besselj(m+1,R1.*beta1))./2.*bessely(m,beta1.*r)).^2.*r;
        C=sqrt(2*pi*integral(g,R1,R1*(1+h1),'AbsTol',0,'RelTol',1e-8));
        vec1(n+1)=d3/(2*d1)/s1*lambda1.^(3/4)*exp(-1i*k*d3*l/2)*R1*h1*integral2(f1,0,1,0,2*pi)/C;     
    end
    vec1(1)=vec1(1)*s_mn_N(k,m,0,R1,h1)/eigenvalue_N(m,0,R1,h1).^(3/4);
    if(mode==1)
        vec1=1i*vec1;
    end
    vec=[vec1];
end
end