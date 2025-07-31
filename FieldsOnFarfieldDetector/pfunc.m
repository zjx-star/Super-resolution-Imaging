function vec = pfunc(c,X,Y,Z,h,ang)
%UNTITLED29 Summary of this function goes here
%   Detailed explanation goes here

vec=0*X;
n=length(X);
X=X+(pi*sqrt(2)-pi/2);
Y=Y-(pi*sqrt(2)-pi/2);
for i=1:n
    z=X(i)+Y(i)*1i;
    a=angle(z)+ang;
    rho=abs(z);
    X(i)=rho*cos(a);
    Y(i)=rho*sin(a);
    if(abs(X(i))>=h/2-h*h/1000 && abs(X(i))<=pi/3+h/2+h*h/1000 && abs(Y(i))<=pi/6+h*h/1000)
        vec(i)=c;
    else
        vec(i)=0;
    end
end

end