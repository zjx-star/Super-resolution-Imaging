function val = G_H_m(beta,r1,r2,order,terms)
% Compute the Fourier coefficients of the fundamental solution
%Input:
%Precondition:
%Output:
%Postcondition:

k=2.*sqrt(r1.*r2)./(r1+r2);
gam=beta.*(r1+r2);
val1=zeros(size(r1,1),size(r2,2));
val2=zeros(size(r1,1),size(r2,2));
for n=[0:terms-1]
    val1=gamma(order-n+1/2).*(gam.^2./4).^(n)./factorial(n).*hypergeom([order+1/2-n,order+1/2],2*order+1,k.^2)+val1;
end
val1=1./factorial(order)./sqrt(pi.*r1.*r2).*(k./2).^(2*order+1).*val1;
for n=[0:terms-1]
    val2=val2+1./gamma(n+order+3/2)./factorial(n).*(-gam.^2/4).^(n).*hypergeom(order+1/2,[n+order+3/2,2.*order+1],gam.^2.*k.^2./4);
end
val2=1i.*sqrt(pi)./factorial(order)./sqrt(r1.*r2).*(k./2).^(2.*order+1).*val2;
val=val1+val2;
end

