function vec = asymp_C_TM_R(m,N,k,k1,R1,R2,h,mode)
%asymp_C_TM is used to calculate the asymptotic values of C_m_TM
%
%Input:
%Precondition:
%Output:
%Postcondition:
vec=zeros(1,N);
if(R1==R2)
    m1=m;
    m=abs(m);
    for(i=1:N)
        f=@(s1,s2)log(1./abs(s2)).*cos(i.*pi.*(s1+s2))./sqrt(2)./pi;
        c1=@(x)-x;
        c2=@(x)1-x;
        vec(i)=(integral2(f,0,1,c1,0)+integral2(f,0,1,0,c2));
    end
    temp=sqrt([1:N]);
    s=s_mn_N(k1,m,0,R1,h);
    vec=-2*m1*k.^2./k1.^2.*(h/R1).^(1/2)*cos(s/2+pi/2*mode)*sqrt(pi).*vec.*temp;
else
end

end
