function k = kasymptotic(m1,m2,eps_ratio,mu_ratio,l,R,h,mode)
%The function is used to compute the asymptotic value of resonances
%
%Input:
%mode: 0 implies even mode while 1 implies odd mode
%Precondition:
%Output:
%Postcondition:

if(m1==0)
    k_0_m2=(2*m2+(1-(-1).^(mode))/2)*pi/l;
    alpha1=3/2/pi+1/pi*integral(@(x)(cos(2.*k_0_m2.*sin(x))-1).*cos(2.*x)./sin(x),0,pi/2)+1/pi*(log(2)-0.5772156649-psi(3/2));
    beta1=integral(@(x)sin(k_0_m2.*2.*sin(x)).*cos(2.*x)./sin(x),0,pi/2)./pi;
    PI=-h*log(h)/pi+alpha1*h+1i*beta1*h-4*(0.003804423620491)*h;
    k=k_0_m2-2*k_0_m2*PI;
else
    m1=abs(m1);
    k_m1_m2=sqrt((m1/R)^2+(m2*pi/l)^2);
    k0=k_m1_m2/sqrt(eps_ratio)/sqrt(mu_ratio);
    alpha=@(m,x)integral(@(th)(cos(2.*x.*sin(th))-1).*cos(2.*m.*th)./sin(th),0,pi/2)+3/2+log(2)-0.5772156649-psi(m+1/2);
    beta=@(m,x)integral(@(th)sin(2.*x.*sin(th)).*cos(2.*(m).*th)./sin(th),0,pi/2);
    %c1=p_constant(20,mu_ratio)
    %c2=p_constant(20,eps_ratio)
    c1=0.003736769907957;
    c2=0.001874856975762;
    Pi=@(x)m1^2/R/pi*(alpha(m1,x*R)+1i*beta(m1,x*R))-x^2*R/pi*((alpha(m1-1,x*R)+alpha(m1+1,x*R))/2+1i*(beta(m1-1,x*R)+beta(m1+1,x*R))/2)-4*m1^2/R*c1+4*x^2*R*c2;
    if(m2==0)
        k=k_m1_m2-(1/mu_ratio*m1^2/R^2-k0^2)*R/pi/l/k_m1_m2*h*log(h)-1/mu_ratio*m1^2/R^2/k_m1_m2/2*h+1/mu_ratio/l/k_m1_m2*Pi(k0)*h;
    else
        k=k_m1_m2-2*(1/mu_ratio*m1^2/R^2-k0^2)*R/pi/l/k_m1_m2*h*log(h)-1/mu_ratio*m1^2/R^2/k_m1_m2/2*h+2/mu_ratio/l/k_m1_m2*Pi(k0)*h;
    end
end
end

