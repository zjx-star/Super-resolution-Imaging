h=0.1;
m=1;
n=0;
mode=0;
Out=[];
format long;
for hs=[h,h/2,h/2.^2,h/2.^3]
    k=kasymptotic(m,n,1,hs,mode);
    k0=real(k);
    D1=D_m(m,1,k,hs,mode);
    D2=D_m(m,1,k0,hs,mode);
    delta_D1=(D1-D2);
    beta_m=1/(4*pi)*(integral(@(x)sin(k0.*2.*sin(x)).*cos(2.*m.*x)./sin(x),0,pi/2));
    beta_m1=1/(4*pi)*(integral(@(x)sin(k0.*2.*sin(x)).*cos(2.*(m+1).*x)./sin(x),0,pi/2));
    beta_m2=1/(4*pi)*(integral(@(x)sin(k0.*2.*sin(x)).*cos(2.*(m-1).*x)./sin(x),0,pi/2));
    delta_D2=4*(k0^2*(beta_m1+beta_m2)/2-m^2*beta_m)*1i*hs;
    delta=delta_D1-delta_D2;
    Out=[Out;delta_D1,delta_D2,((delta(1:end)).*((delta(1:end)).')').^(1/2)];
end
display(Out);