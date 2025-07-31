function vec = g21(k,l,tar,X,Y,Z)
%g11 Summary of this function goes here
%Generate the vector to compute the g11 field at tar 

x=tar(1);
y=tar(2);
z=tar(3);
dis=sqrt((X-x).^2+(Y-y).^2+(Z-z).^2);
vec=(-k^2*exp(1i*k*dis)./(4*pi*dis.^3)-3i*k*exp(1i*k*dis)./(4*pi*dis.^4)+3*exp(1i*k*dis)./(4*pi*dis.^5)).*(Y-y).*(X-x);

dis=sqrt((X-x).^2+(Y-y).^2+(Z+z-l).^2);
vec=vec-((-k^2*exp(1i*k*dis)./(4*pi*dis.^3)-3i*k*exp(1i*k*dis)./(4*pi*dis.^4)+3*exp(1i*k*dis)./(4*pi*dis.^5)).*(Y-y).*(X-x));
end