function vec = asymp_R_TE_R(m,N,k,k1,R1,R2,h,mode)
%asymp_R_TM is used to calculate the asymptotic values of R_m_TE
%
%Input:
%Precondition:
%Output:
%Postcondition:
vec=zeros(1,N);
m=abs(m);
if(R1==R2)
    for(i=1:N)
        f=@(s1,s2)log(1./abs(s2)).*cos(i.*pi.*(s1+s2))./sqrt(2)./pi;
        c1=@(x)-x;
        c2=@(x)1-x;
        vec(i)=-2i*(R1*h).^(1/2).*(integral2(f,0,1,c1,0)+integral2(f,0,1,0,c2));
    end
    temp=sqrt([1:N]);
    vec=sqrt(pi).*vec.*temp;
else

end
end

