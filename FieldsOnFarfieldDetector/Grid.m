function [X,Y,Z,w] = Grid(point1,point2,Nx,Ny,Nz)
%Generate the grids and weights for computing far field 

x1=point1(1);
y1=point1(2);
z1=point1(3);
x2=point2(1);
y2=point2(2);
z2=point2(3);

X=[];
Y=[];
Z=[];
wX=[];
wY=[];
wZ=[];

[x,w]=lgwt(4,-1,1);
x=x(end:-1:1);
for i=1:Nx
    a=x1+(i-1)*(x2-x1)/Nx;
    b=x1+i*(x2-x1)/Nx;
    X=[X,(b-a)/2*x.'+(a+b)/2];
    wX=[wX,(b-a)/2*w.'];
end

for i=1:Ny
    a=y1+(i-1)*(y2-y1)/Ny;
    b=y1+i*(y2-y1)/Ny;
    Y=[Y,(b-a)/2*x.'+(a+b)/2];
    wY=[wY,(b-a)/2*w.'];
end

[x,w]=lgwt(2,-1,1);
x=x(end:-1:1);
for i=1:Nz
    a=z1+(i-1)*(z2-z1)/Nz;
    b=z1+i*(z2-z1)/Nz;
    Z=[Z,(b-a)/2*x.'+(a+b)/2];
    wZ=[wZ,(b-a)/2*w.'];
end

[X,Y,Z]=ndgrid(X,Y,Z);  
[wX,wY,wZ]=ndgrid(wX,wY,wZ);

X=reshape(X,[1,Nx*Ny*Nz*4^2*2]);
Y=reshape(Y,[1,Nx*Ny*Nz*4^2*2]);
Z=reshape(Z,[1,Nx*Ny*Nz*4^2*2]);
wX=reshape(wX,[1,Nx*Ny*Nz*4^2*2]);
wY=reshape(wY,[1,Nx*Ny*Nz*4^2*2]);
wZ=reshape(wZ,[1,Nx*Ny*Nz*4^2*2]);
w=wX.*wY.*wZ;
end