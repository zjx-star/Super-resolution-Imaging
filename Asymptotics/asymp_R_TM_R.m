function vec = asymp_R_TM_R(m,N,k,k1,R1,R2,h,mode)
%asymp_R_TM is used to calculate the asymptotic values of R_m_TM
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
    lambda=eigenvalue_N(m,0,R1,h);
    temp=sqrt([1:N]);
    vec=-2.*k.^2.*m1./lambda*h.^(1/2)/R1.^(1/2).*sqrt(pi).*vec.*temp;
else

end
end

