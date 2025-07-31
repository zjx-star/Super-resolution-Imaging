function [X,Y,Z,w] = Grid2(point1,point2,Nx,Ny,Nz)
%Generate the grids and weights for computing far field 

x1=point1(1);
y1=point1(2);
z1=point1(3);
x2=point2(1);
y2=point2(2);
z2=point2(3);


X=[0:Nx-1]*(x2-x1)/(Nx-1)+x1;
Y=[0:Ny-1]*(y2-y1)/(Ny-1)+y1;
Z=[0:Nz-1]*(z2-z1)/(Nz-1)+z1;
wX=0*X+(x2-x1)/(Nx-1);
wX([1,Nx])=wX([1,Nx])/2;
wY=0*Y+(y2-y1)/(Ny-1);
wY([1,Ny])=wY([1,Ny])/2;
wZ=0*Z+(z2-z1)/(Nz-1);
wZ([1,Nz])=wZ([1,Nz])/2;

[X,Y,Z]=ndgrid(X,Y,Z);  
[wX,wY,wZ]=ndgrid(wX,wY,wZ);

X=reshape(X,[1,Nx*Ny*Nz]);
Y=reshape(Y,[1,Nx*Ny*Nz]);
Z=reshape(Z,[1,Nx*Ny*Nz]);
wX=reshape(wX,[1,Nx*Ny*Nz]);
wY=reshape(wY,[1,Nx*Ny*Nz]);
wZ=reshape(wZ,[1,Nx*Ny*Nz]);
w=wX.*wY.*wZ;
end