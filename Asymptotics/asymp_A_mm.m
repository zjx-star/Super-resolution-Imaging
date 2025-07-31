function val = asymp_A_mm(m,k,k1,h,mode)
%Calculate the asymptotic value of A_mm
%
%Input:
%mode: 0 implies even mode while 1 implies odd mode
%Precondition:
%Output:
%Postcondition:

if(m==0)
    val=-log(h)+3/2+log(2)-0.5772156649-psi(3/2)+integral(@(x)(cos(2.*k.*sin(x))-1).*cos(2.*x)./sin(x),0,pi/2)+integral(@(x)sin(2.*k.*sin(x)).*cos(2.*x)./sin(x),0,pi/2)*1i;
    val=-val*k*cos(k/2+pi/2*mode)*h/pi;
else
    m=abs(m);
    Gamma1=-log(h)+3/2+log(2)-0.5772156649-psi(m+1/2)+integral(@(x)(cos(2.*k.*sin(x))-1).*cos(2.*m.*x)./sin(x),0,pi/2)+integral(@(x)sin(2.*k.*sin(x)).*cos(2.*m.*x)./sin(x),0,pi/2)*1i;
    Gamma2=-log(h)+3/2+log(2)-0.5772156649-psi(m+1+1/2)+integral(@(x)(cos(2.*k.*sin(x))-1).*cos(2.*(m+1).*x)./sin(x),0,pi/2)+integral(@(x)sin(2.*k.*sin(x)).*cos(2.*(m+1).*x)./sin(x),0,pi/2)*1i;
    Gamma3=-log(h)+3/2+log(2)-0.5772156649-psi(m-1+1/2)+integral(@(x)(cos(2.*k.*sin(x))-1).*cos(2.*(m-1).*x)./sin(x),0,pi/2)+integral(@(x)sin(2.*k.*sin(x)).*cos(2.*(m-1).*x)./sin(x),0,pi/2)*1i;
    s=s_mn_N(k,m,0,1,h);
    lambda=eigenvalue_N(m,0,1,h);
    val=(2*cos(s/2+pi/2*mode)*lambda*Gamma1-m.^2*k.^2*cos(s/2+pi/2*mode)/lambda*(Gamma2+Gamma3))*h/(2*pi);
end
end

