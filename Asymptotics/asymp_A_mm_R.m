function val = asymp_A_mm_R(m,k,k1,R1,R2,h,mode)
%Calculate the asymptotic value of A_mm
%
%Input:
%mode: 0 implies even mode while 1 implies odd mode
%Precondition:
%Output:
%Postcondition:

if(R1==R2)
    m=abs(m);
    s=s_mn_N(k1,m,0,R1,h);
    lambda=eigenvalue_N(m,0,R1,h);


    Gamma1=3/2+log(2)-0.5772156649-psi(m+1/2)+integral(@(x)(cos(2.*k.*R1.*sin(x))-1).*cos(2.*m.*x)./sin(x),0,pi/2)+integral(@(x)sin(2.*k.*R1.*sin(x)).*cos(2.*m.*x)./sin(x),0,pi/2)*1i;
    Gamma2=3/2+log(2)-0.5772156649-psi(m+1+1/2)+integral(@(x)(cos(2.*k.*R1.*sin(x))-1).*cos(2.*(m+1).*x)./sin(x),0,pi/2)+integral(@(x)sin(2.*k.*R1.*sin(x)).*cos(2.*(m+1).*x)./sin(x),0,pi/2)*1i;
    Gamma3=3/2+log(2)-0.5772156649-psi(m-1+1/2)+integral(@(x)(cos(2.*k.*R1.*sin(x))-1).*cos(2.*(m-1).*x)./sin(x),0,pi/2)+integral(@(x)sin(2.*k.*R1.*sin(x)).*cos(2.*(m-1).*x)./sin(x),0,pi/2)*1i;
    
    val=-(lambda-m^2*k^2/lambda/R1^2)*cos(s/2)*R1*h*log(h)/pi+cos(s/2)*lambda*R1*h/pi*Gamma1-m^2*k^2*cos(s/2)*R1*h/(2*lambda*pi)*(Gamma2+Gamma3)/R1^2;
else
    s=s_mn_N(k1,m,0,R2,h);
    Gamma1=integral(@(x)(exp(1i.*k.*sqrt(R1.^2+R2.^2-2.*R1.*R2.*cos(x)))./sqrt(R1.^2+R2.^2-2.*R1.*R2.*cos(x)).*cos(m.*x)),0,2*pi);
    Gamma2=integral(@(x)(exp(1i.*k.*sqrt(R1.^2+R2.^2-2.*R1.*R2.*cos(x)))./sqrt(R1.^2+R2.^2-2.*R1.*R2.*cos(x)).*cos((m+1).*x)),0,2*pi);
    Gamma3=integral(@(x)(exp(1i.*k.*sqrt(R1.^2+R2.^2-2.*R1.*R2.*cos(x)))./sqrt(R1.^2+R2.^2-2.*R1.*R2.*cos(x)).*cos((m-1).*x)),0,2*pi);
    val=-2*cos(s/2)*(-m^2*h*R1/(4*pi*R2)*Gamma1+R1^2*k^2*h/(8*pi)*(Gamma2+Gamma3));
end
end

