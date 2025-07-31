function val = p_constant2(N)
%p_constant is used to calculate the parameter p(1+4*P)^(-1)p
%
%Input:
%Precondition:
%Output:
%Postcondition:
p=zeros(N,1);

f=@(s1,s2)log(1./sqrt(abs(s2)))./pi;
c1=@(x)-x;
c2=@(x)1-x;
a=2*(integral2(f,0,1,c1,0)+integral2(f,0,1,0,c2));
a=0;

for(i=1:N)
    f=@(s1,s2)log(1./abs(s2)).*cos(i.*pi.*(s1+s2))./sqrt(2)./pi;
    c1=@(x)-x;
    c2=@(x)1-x;
    p(i)=sqrt(i*pi)*(integral2(f,0,1,c1,0)+integral2(f,0,1,0,c2));
end

P=zeros(N,N);
for(i=1:N)
    for(j=1:N)
       f=@(s1,s2)log(1./abs(s2)).*cos(j.*pi.*s1).*cos(i.*pi.*(s1+s2))./2./pi;
       c1=@(x)-x;
       c2=@(x)1-x;
       P(i,j)=sqrt(i*j)*pi*(integral2(f,0,1,c1,0)+integral2(f,0,1,0,c2)); 
    end
end

P=2*P;
P(1:N+1:N*N)=P(1:N+1:N*N)+1;

A=zeros(N+1,N+1);
A(1,1)=a;
A(1,2:N+1)=2*p.';
A(2:N+1,1)=2*p;
A(2:N+1,2:N+1)=P;
eig(A)
phi=zeros(N+1,1);
phi(1)=1;
val=(A(1,1)-1/((phi.')*(A\phi)))/4;
X=A\phi;
X'
[B,C]=eig(A);
B
C
%val = -1/((phi.')*(A\phi))/4;



end

