kinc=1;
h=0.1;
l=1;
x=[-pi:0.01:pi];
y=[-pi:0.01:pi];
len=length(x);
[X,Y]=meshgrid(x,y);
X=X(:);
Y=Y(:);
Z=0*X+l;


Psi=@(x1,x2,x3,y1,y2)exp(1i.*kinc.*sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2))./4/pi./sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2);
Partial_x1_Psi=@(x1,x2,x3,y1,y2)(1i.*kinc.*exp(1i.*kinc.*sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2))./4/pi./sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2)-...
    exp(1i.*kinc.*sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2))./4/pi./(x3.^2+(x1-y1).^2+(x2-y2).^2)).*(x1-y1)./sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2);
Partial_x2_Psi=@(x1,x2,x3,y1,y2)(1i.*kinc.*exp(1i.*kinc.*sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2))./4/pi./sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2)-...
    exp(1i.*kinc.*sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2))./4/pi./(x3.^2+(x1-y1).^2+(x2-y2).^2)).*(x2-y2)./sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2);
Partial_x3_Psi=@(x1,x2,x3,y1,y2)(1i.*kinc.*exp(1i.*kinc.*sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2))./4/pi./sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2)-...
    exp(1i.*kinc.*sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2))./4/pi./(x3.^2+(x1-y1).^2+(x2-y2).^2)).*(x3)./sqrt(x3.^2+(x1-y1).^2+(x2-y2).^2);


integrand_Psi = @(x1,x2,x3,s,theta)h.*Psi(x1,x2,x3,(1+s.*h).*cos(theta),(1+s.*h).*sin(theta)).*(1+s.*h);
Psi = AnnularFarField(integrand_Psi,X,Y,Z);
X=reshape(X,len,len);
Y=reshape(Y,len,len);
Z1=abs(reshape(Psi,len,len));
mesh(X,Y,Z1);