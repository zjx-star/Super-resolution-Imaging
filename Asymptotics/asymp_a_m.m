function val = asymp_a_m(m,k,l,r,h,d1,mode)
%Calculate the asymptotic value of am
%
%Input:
%Precondition:
%Output:
%Postcondition:

if(nargin<6 || d1==0)
    lambda=eigenvalue_N(m,0,h);
    val=-k*m.^(3/2)*exp(-1i*k*l/2)*sqrt(2*pi*h)/(4i);
    if(mode==1)
        val=1i*val;
    end
else
    d3=sqrt(1-d1.^2);
    lambda=eigenvalue_N(m,0,r,h);
    val=(1i).^(m)/2*r*besselj(m,k*d1*r)/d1*d3*exp(-1i*k*d3*l/2)*sqrt(2*pi*h);
    if(mode==1)
        val=1i*val;
    end
end
end

