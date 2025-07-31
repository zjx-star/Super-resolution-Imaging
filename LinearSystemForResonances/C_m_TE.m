function vec = C_m_TE(m,N,l,k,k1,R1,R2,h1,h2,mode)
%Compute C_m^TE
%
%Input:
%mode: 0 imples even mode while 1 implies odd mode
%Precondition:
%Output:
%Postcondition:

a=R1*h1/(R2*h2);
vec=zeros(1,N);
if(m==0)
    %parfor n=1:N
    %    beta=accu_beta_N(0,n,h)
    %    lambda=eigenvalue_N(0,n,h)
    %    s=s_mn_N(k,0,n,h)
    %    f=@(s1,s2,theta)exp(1i.*k.*sqrt(h.^2.*s2.^2+4.*(1+s1.*h).*(1+(s2+s1).*h).*sin(theta./2).^2))./(4.*pi.*sqrt(h.^2.*s2.^2+4.*(1+s1.*h).*(1+(s2+s1).*h).*sin(theta./2).^2)).*...
    %        ((bessely(m-1,beta)-bessely(m+1,beta))./2.*(besselj(m-1,beta.*(1+(s2+s1).*h))-besselj(m+1,beta.*(1+(s2+s1).*h)))./2-(besselj(m-1,beta)-besselj(m+1,beta))./2.*(bessely(m-1,beta.*(1+(s2+s1).*h))-bessely(m+1,beta.*(1+(s2+s1).*h)))./2).*...
    %        sin(theta).*(1+(s1+s2).*h);
    %    g=@(r)((bessely(m-1,beta)-bessely(m+1,beta))./2.*besselj(m,beta.*r)-(besselj(m-1,beta)-besselj(m+1,beta))./2.*bessely(m,beta.*r)).^2;
    %    C=sqrt(2*pi*integral(g,1,1+h,'AbsTol',1e-8,'RelTol',1e-8));
    %    c1=@(x)-x;
    %    c2=@(x)1-x;
    %    vec(n)=-2*cos(k*l/2)./lambda.^(1/4)./s.*...
    %        2*pi*k^2*h.^2.*beta.*2*integral3(f,0,1,c1,c2,0,pi)/...
    %        C;        
    %end
else
    m1=m;
    m=abs(m);
    for n = 1:N
        beta1=accu_beta_N(m,n,R1,h1);
        beta2=accu_beta_N(m,0,R2,h2);
        lambda1=eigenvalue_N(m,n,R1,h1);
        lambda2=eigenvalue_N(m,0,R2,h2);
        s1=s_mn_N(k1,m,n,R1,h1);
        s2=s_mn_N(k1,m,0,R2,h2);
        f1=@(s1,s2,theta)exp(1i.*k.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2))./(4.*pi.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2)).*...
            ((bessely(m-1,R1*beta1)-bessely(m+1,R1*beta1))./2.*besselj(m,beta1.*R1.*(1+s1.*h1))-(besselj(m-1,R1*beta1)-besselj(m+1,R1*beta1))./2.*bessely(m,beta1.*R1.*(1+s1.*h1))).*...
            ((bessely(m-1,R2*beta2)-bessely(m+1,R2*beta2))./2.*besselj(m,beta2.*R2.*(1+(s2+a*s1).*h2))-(besselj(m-1,R2*beta2)-besselj(m+1,R2*beta2))./2.*bessely(m,beta2.*R2.*(1+(s2+a*s1).*h2))).*...
            cos(m1.*theta).*R1.*(1+s1.*h1).*R2.*(1+(a*s1+s2).*h2);
        f2=@(s1,s2,theta)exp(1i.*k.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2))./(4.*pi.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2)).*...
            ((bessely(m-1,R1*beta1)-bessely(m+1,R1*beta1))./2.*(besselj(m-1,beta1.*R1.*(1+s1.*h1))-besselj(m+1,beta1.*R1.*(1+s1.*h1)))./2-(besselj(m-1,R1.*beta1)-besselj(m+1,R1.*beta1))./2.*(bessely(m-1,beta1.*R1.*(1+s1.*h1))-bessely(m+1,beta1.*R1.*(1+s1.*h1)))./2).*...
            ((bessely(m-1,R2*beta2)-bessely(m+1,R2*beta2))./2.*(besselj(m-1,beta2.*R2.*(1+(s2+a.*s1).*h2))-besselj(m+1,beta2.*R2.*(1+(s2+a.*s1).*h2)))./2-(besselj(m-1,R2.*beta2)-besselj(m+1,R2.*beta2))./2.*(bessely(m-1,beta2.*R2.*(1+(s2+a.*s1).*h2))-bessely(m+1,beta2.*R2.*(1+(s2+a.*s1).*h2)))./2).*...
            cos(m1.*theta).*cos(theta).*R1.*(1+s1.*h1).*R2.*(1+(a.*s1+s2).*h2);
        f3=@(s1,s2,theta)exp(1i.*k.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2))./(4.*pi.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2)).*...
            ((bessely(m-1,R1*beta1)-bessely(m+1,R1*beta1))./2.*besselj(m,beta1.*R1.*(1+s1.*h1))-(besselj(m-1,R1.*beta1)-besselj(m+1,R1.*beta1))./2.*bessely(m,beta1.*R1.*(1+s1.*h1))).*...
            ((bessely(m-1,R2*beta2)-bessely(m+1,R2*beta2))./2.*(besselj(m-1,beta2.*R2.*(1+(s2+a.*s1).*h2))-besselj(m+1,beta2.*R2.*(1+(s2+a.*s1).*h2)))./2-(besselj(m-1,R2.*beta2)-besselj(m+1,R2.*beta2))./2.*(bessely(m-1,beta2.*R2.*(1+(s2+a.*s1).*h2))-bessely(m+1,beta2.*R2.*(1+(s2+a.*s1).*h2)))./2).*...
            -1i.*sin(m1.*theta).*sin(theta).*R2.*(1+(a.*s1+s2).*h2);
        f4=@(s1,s2,theta)exp(1i.*k.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2))./(4.*pi.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2)).*...
            ((bessely(m-1,R1*beta1)-bessely(m+1,R1*beta1))./2.*(besselj(m-1,beta1.*R1.*(1+s1.*h1))-besselj(m+1,beta1.*(1+s1.*h1)))./2-(besselj(m-1,beta1)-besselj(m+1,beta1))./2.*(bessely(m-1,beta1.*(1+s1.*h1))-bessely(m+1,beta1.*(1+s1.*h1)))./2).*...
            ((bessely(m-1,R2*beta2)-bessely(m+1,R2*beta2))./2.*besselj(m,beta2.*R2.*(1+(s2+a.*s1).*h2))-(besselj(m-1,R2.*beta2)-besselj(m+1,R2.*beta2))./2.*bessely(m,beta2.*R2.*(1+(s2+a.*s1).*h2))).*...
            -1i.*sin(m1.*theta).*sin(theta).*R1.*(1+s1.*h1);
        f5=@(s1,s2,theta)exp(1i.*k.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2))./(4.*pi.*sqrt((R2-R1+R2.*h2.*s2).^2+4.*R1.*R2.*(1+h1.*s1).*(1+h2.*(s2+a.*s1)).*sin(theta/2).^2)).*...
            ((bessely(m-1,R1*beta1)-bessely(m+1,R1*beta1))./2.*besselj(m,beta1.*R1.*(1+s1.*h1))-(besselj(m-1,R1*beta1)-besselj(m+1,R1*beta1))./2.*bessely(m,beta1.*R1.*(1+s1.*h1))).*...
            ((bessely(m-1,R2*beta2)-bessely(m+1,R2*beta2))./2.*besselj(m,beta2.*R2.*(1+(s2+a.*s1).*h2))-(besselj(m-1,R2*beta2)-besselj(m+1,R2*beta2))./2.*bessely(m,beta2.*R2.*(1+(s2+a.*s1).*h2))).*...
            cos(m1.*theta).*cos(theta);
        g1=@(r)((bessely(m-1,R1*beta1)-bessely(m+1,R1*beta1))./2.*besselj(m,beta1.*r)-(besselj(m-1,R1*beta1)-besselj(m+1,R1*beta1))./2.*bessely(m,beta1.*r)).^2.*r;
        g2=@(r)((bessely(m-1,R2*beta2)-bessely(m+1,R2*beta2))./2.*besselj(m,beta2.*r)-(besselj(m-1,R2*beta2)-besselj(m+1,R2*beta2))./2.*bessely(m,beta2.*r)).^2.*r;
        C1=sqrt(2*pi*integral(g1,R1,R1*(1+h1),'AbsTol',1e-8,'RelTol',1e-8));
        C2=sqrt(2*pi*integral(g2,R2,R2*(1+h2),'AbsTol',1e-8,'RelTol',1e-8));
        c1=@(x)-a*x;
        c2=@(x)1-a*x;
        vec(n)=-2*cos(s2*l/2+pi/2*mode)./s1/lambda1.^(1/4)/C1./C2.*...
            2.*pi.*R1*h1*R2*h2.*2.*(-lambda1.*lambda2.*integral3(f1,0,1,c1,c2,0,pi)+k.^2.*beta1.*beta2.*integral3(f2,0,1,c1,c2,0,pi)-k.^2.*beta2.*(-1i.*m1).*integral3(f3,0,1,c1,c2,0,pi)+...
            k.^2.*beta1.*(1i.*m1).*integral3(f4,0,1,c1,c2,0,pi)+k.^2.*m1.^2.*integral3(f5,0,1,c1,c2,0,pi));
    end
end

end

