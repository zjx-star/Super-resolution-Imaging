function val = accu_beta_N(m,n,R,h,delta,epsilon,maxi)
%Compute betta_N defined in eq. (3.11)
%
%Input:
%m: angular momentum of the eigenvalue defined in eq. (3.11)
%n: the longitudinal index of the eigenvalue defined in eq. (3.11)
%R,h: the eigenvalue is in an annular aperture of radius R and width R*h
%delta,epsilon,maxi: stop criteria for computing beta_N by Newton's method
%
%Output:
%val: value of betta_N defiend in eq. (3.11)
%
if(nargin<5)
    delta=1e-10;
    epsilon=1e-10;
    maxi=100;
end

m=abs(m);

f=@(x)(bessely(m-1,R*x)-bessely(m+1,R*x))/2*(besselj(m-1,(1+h)*R*x)-besselj(m+1,(1+h)*R*x))/2-(besselj(m-1,R*x)-besselj(m+1,R*x))/2*(bessely(m-1,(1+h)*R*x)-bessely(m+1,(1+h)*R*x))/2;
df=@(x)1/4*R*(bessely(m-2,R*x)-2*bessely(m,R*x)+bessely(m+2,R*x))*1/2*(besselj(m-1,(1+h)*R*x)-besselj(m+1,(1+h)*R*x)) +...
    1/2*(1+h)*(bessely(m-1,R*x)-bessely(m+1,R*x))*1/4*R*(besselj(m-2,(1+h)*R*x)-2*besselj(m,(1+h)*R*x)+besselj(m+2,(1+h)*R*x)) -...
    1/4*R*(besselj(m-2,R*x)-2*besselj(m,R*x)+besselj(m+2,R*x))*1/2*(bessely(m-1,(1+h)*R*x)-bessely(m+1,(1+h)*R*x)) -...
    1/2*(1+h)*(besselj(m-1,R*x)-besselj(m+1,R*x))*1/4*R*(bessely(m-2,(1+h)*R*x)-2*bessely(m,(1+h)*R*x)+bessely(m+2,(1+h)*R*x));
x_asymp=asymp_beta_N(m,n,R,h);
x_accu=newton(f,df,x_asymp,delta,epsilon,maxi);
val=x_accu;
end

