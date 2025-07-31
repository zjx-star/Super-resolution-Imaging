function [p,y,err] = muller(f,p0,p1,p2,delta,epsilon,max1)
%Input  - f is the object function input as a string 'f'
%       - p0,p1, and p2 are the initial approximations
%       - delta is the tolerance for p0,p1, and p2
%       - epsilon the tolerance for the function values y
%       - max1 is the maximum number of iterations
%Output - p is the Muller approximation to the zero of f
%       - y is the function value y=f(p) 
%       - err is the error in the approximation of p.

%Initialize the matrices P and Y
P=[p0 p1 p2];
Y=[f(P(1)) f(P(2)) f(P(3))];

%Calculate a and b in formula (15)
for k=1:max1
    h0=P(1)-P(3);h1=P(2)-P(3);e0=Y(1)-Y(3);e1=Y(2)-Y(3);c=Y(3);
    denom=h1*h0^2-h0*h1^2;
    a=(e0*h1-e1*h0)/denom;
    b=(e1*h0^2-e0*h1^2)/denom;
    
    %Observe the variation of P and Y
    %disp(P);
    %disp(Y);
    
    %Suppress any complex roots 
    %if b^2-4*a*c > 0
    %    disc=sqrt(b^2-4*a*c);
    %else
    %    disc=0;
    %end
    
    disc=sqrt(b^2-4*a*c);
    
    %Find the smallest root of (17)
    if (b<0) 
        disc=-disc;
    end
    z=-2*c/(b+disc);
    p=P(3)+z;
    
    if(imag(p)>0)
       z=-2*c/(b-disc);
       p=P(3)+z;
    end
    
    %Sort the entries of P to find the two closest to p
    if abs(p-P(2))<abs(p-P(1))
        Q=[P(2) P(1) P(3)];
        P=Q;
        Y=[f(P(1)) f(P(2)) f(P(3))];
    end
    if abs(p-P(3))<abs(p-P(2))
        R=[P(1) P(3) P(2)];
        P=R;
        Y=[f(P(1)) f(P(2)) f(P(3))];
    end
   
    %Replace the entry of P that was farthest from p with p
    P(3)=p;
    Y(3)=feval(f,P(3));
    y=Y(3);
    
    
    %Determine stopping criteria
    err=abs(z);
    relerr=err/(abs(p)+delta);
    if(err<delta) || (relerr<delta) || (abs(y)<epsilon)
        break;
    end
end

