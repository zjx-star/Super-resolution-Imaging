function mat = asymp_B_TE_TE(m,N,k,l,R1,R2,h,mode)
%asymp_B_TE_TE is used to calculate the asymptotic values of asymp_B_TM_TM
%
%Input:
%Precondition:
%Output:
%Postcondition:
mat=zeros(N,N);
if(R1==R2)
    m=abs(m);
    for(i=1:N)
        for(j=1:N)
           f=@(s1,s2)log(1./abs(s2)).*cos(j.*pi.*s1).*cos(i.*pi.*(s1+s2))./2./pi;
           c1=@(x)-x;
           c2=@(x)1-x;
           mat(i,j)=-4.*(i)^(1/2).*(j).^(1/2).*pi.*(integral2(f,0,1,c1,0)+integral2(f,0,1,0,c2)); 
        end
    end   
else
end

end

