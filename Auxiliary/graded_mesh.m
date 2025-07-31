function [w,dw] = graded_mesh(p1,p2,a,b,n)
%This function is used to construct the graded mesh
%   
p=max(p1,p2);
zeta=[0:n-1]*2*pi/n;
ksi=(2*zeta-(2*pi))/(2*pi);
w1=(1/2-1/p)*ksi.^3+ksi./p+1/2;
w2=1-w1;
w=(b.*w1.^p1+a.*w2.^p2)./(w1.^p1+w2.^p2);
dw1=(3*(1/2-1/p).*ksi.^2+1/p).*2./(2.*pi);
dw2=-dw1;
dw=(p1.*b.*w1.^(p1-1).*dw1+p2.*a.*w2.^(p2-1).*dw2)./(w1.^p1+w2.^p2)-(b.*w1.^p1+a.*w2.^p2).*(p1.*w1.^(p1-1).*dw1+p2.*w2.^(p2-1).*dw2)./(w1.^p1+w2.^p2).^2;
end