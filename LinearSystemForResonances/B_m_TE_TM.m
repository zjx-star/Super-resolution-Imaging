function mat = B_m_TE_TM(m,N,l,k,k1,R1,R2,h1,h2,mode)
%Construct B_m_TE_TE
%   
%Input:
%mode: 0 imples even mode while 1 imples odd mode
%Precondition:
%Output
%Postcondition

a=R1*h1/(R2*h2);
mat=zeros(N,N);
if(m==0)
    
else
    m1=m;
    m=abs(m);
    for(i=0:N*N-1)
        n1=floor(i/N)+1;
        n2=mod(i,N)+1;
        beta1=accu_beta_N(m,n1,R1,h1);
        beta2=accu_beta_D(m,n2,R2,h2);
        s1=s_mn_N(k1,m,n1,R1,h1);
        s2=s_ij_D(k1,m,n2,R2,h2);
        lambda1=eigenvalue_N(m,n1,R1,h1);
        lambda2=eigenvalue_D(m,n2,R2,h2);
        f1=@(s1,s2,theta)exp(1i.*k.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2))./(4.*pi.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2)).*...
            ((bessely(m-1,R1*beta1)-bessely(m+1,R1*beta1))./2.*(besselj(m-1,beta1.*R1.*(1+s1.*h1))-besselj(m+1,beta1.*R1.*(1+s1.*h1)))./2-(besselj(m-1,R1.*beta1)-besselj(m+1,R1.*beta1))./2.*(bessely(m-1,beta1.*R1.*(1+s1.*h1))-bessely(m+1,beta1.*R1.*(1+s1.*h1)))./2).*...
            (bessely(m,R2.*beta2).*(besselj(m-1,R2.*beta2.*(1+(s2+a*s1).*h2))-besselj(m+1,beta2.*(1+(s2+a*s1).*h2)))./2-besselj(m,R2.*beta2).*(bessely(m-1,beta2.*R2.*(1+(s2+a*s1).*h2))-bessely(m+1,beta2.*R2.*(1+(s2+a*s1).*h2)))./2).*...
            -1i.*sin(m1.*theta).*sin(theta).*R1.*(1+s1.*h1).*R2.*(1+(a*s1+s2).*h2);
        f2=@(s1,s2,theta)exp(1i.*k.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2))./(4.*pi.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2)).*...
            ((bessely(m-1,R1*beta1)-bessely(m+1,R1*beta1))./2.*besselj(m,beta1.*R1.*(1+s1.*h1))-(besselj(m-1,R1.*beta1)-besselj(m+1,R1.*beta1))./2.*bessely(m,beta1.*R1.*(1+s1.*h1))).*...
            (bessely(m,R2.*beta2).*(besselj(m-1,beta2.*R2.*(1+(s2+a*s1).*h2))-besselj(m+1,beta2.*R2.*(1+(s2+a*s1).*h2)))./2-besselj(m,R2.*beta2).*(bessely(m-1,beta2.*R2.*(1+(s2+a.*s1).*h2))-bessely(m+1,beta2.*R2.*(1+(s2+a.*s1).*h2)))./2).*...
            cos(m1.*theta).*cos(theta).*R2.*(1+(a*s1+s2).*h2);
        f3=@(s1,s2,theta)exp(1i.*k.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2))./(4.*pi.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2)).*...
            ((bessely(m-1,R1.*beta1)-bessely(m+1,R1.*beta1))./2.*(besselj(m-1,beta1.*R1.*(1+s1.*h1))-besselj(m+1,beta1.*R1.*(1+s1.*h1)))./2-(besselj(m-1,R1.*beta1)-besselj(m+1,R1.*beta1))./2.*(bessely(m-1,beta1.*R1.*(1+s1.*h1))-bessely(m+1,beta1.*R1.*(1+s1.*h1)))./2).*...
            (bessely(m,R2.*beta2).*besselj(m,beta2.*R2.*(1+(s2+a*s1).*h2))-besselj(m,R2.*beta2).*bessely(m,beta2.*R2.*(1+(s2+a*s1).*h2))).*...
            cos(m1.*theta).*cos(theta).*R1.*(1+(s1).*h1);
        f4=@(s1,s2,theta)exp(1i.*k.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2))./(4.*pi.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2)).*...
            ((bessely(m-1,R1.*beta1)-bessely(m+1,R1.*beta1))./2.*besselj(m,beta1.*R1.*(1+s1.*h1))-(besselj(m-1,R1.*beta1)-besselj(m+1,R1.*beta1))./2.*bessely(m,beta1.*R1.*(1+s1.*h1))).*...
            (bessely(m,R2.*beta2).*besselj(m,beta2.*R2.*(1+(s2+a*s1).*h2))-besselj(m,R2.*beta2).*bessely(m,beta2.*R2.*(1+(s2+a.*s1).*h2))).*...
            -1i.*sin(m1.*theta).*sin(theta);
        g1=@(r)((bessely(m-1,R1*beta1)-bessely(m+1,R1*beta1))./2.*besselj(m,beta1.*r)-(besselj(m-1,R1*beta1)-besselj(m+1,R1*beta1))./2.*bessely(m,beta1.*r)).^2.*r;        
        g2=@(r)(bessely(m,R2*beta2).*besselj(m,beta2.*r)-besselj(m,R2*beta2).*bessely(m,beta2.*r)).^2.*r;
        C1=sqrt(2*pi*integral(g1,R1,R1*(1+h1),'AbsTol',1e-8,'RelTol',1e-8));
        C2=sqrt(2*pi*integral(g2,R2,R2*(1+h2),'AbsTol',1e-8,'RelTol',1e-8));
        c1=@(x)-a*x;
        c2=@(x)1-a*x;
        mat(i+1)=s1^(-1)*-2*cot(s2*l/2+pi/2*mode)/lambda2^(1/4)/lambda1^(1/4)./C1./C2.*...
            2*pi*R1*h1*R2*h2.*2*(k.^2.*beta1.*beta2.*integral3(f1,0,1,c1,c2,0,pi)+k.^2.*(-1i.*m1).*beta2.*integral3(f2,0,1,c1,c2,0,pi)+...
            -k.^2.*(1i.*m1).*beta1.*integral3(f3,0,1,c1,c2,0,pi)+k.^2.*(m1.^2).*integral3(f4,0,1,c1,c2,0,pi));
    end
end


end

