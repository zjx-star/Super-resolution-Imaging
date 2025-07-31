function val = p_constant(N,l)
%p_constant is used to calculate the parameter p(1+4*P)^(-1)p
%
%Input:
%Precondition:
%Output:
%Postcondition:
p=zeros(N,1);

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
P=P*4;
P(1:N+1:N*N)=P(1:N+1:N*N)+l;
val=(p')*(P\p);
end

