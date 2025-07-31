function field = AnnularNearField(f,X,Y,Z)
% AnnularFarField computes the far field by Green's representation 
% theorem 
% 
% Input:
% f: fucntion handle
% X,Y,Z: positions

N=48;
p=6;
x=[0.994700467495825,0.972287511536616,0.932815601193916,0.877702204177502,0.808938122201322,...
0.729008388828614,0.640801775389630,0.547506254918819,0.452493745081181,0.359198224610371,...
0.270991611171386,0.191061877798678,0.122297795822499,0.0671843988060841,0.0277124884633837,0.00529953250417503].';
[w,dw]=graded_mesh(p,p,0,2*pi,N);
[x,w]=ndgrid(x,w);
len=length(X);
[temp,dw]=ndgrid(ones(1,16),dw);
tensor=zeros(16,N,len);
Ang=angle(X+1i*Y);
for(i=1:len)
    tensor(:,:,i)=f(X(i),Y(i),Z(i),x,w+Ang(i)).*dw;
end
field=quad2d_legendre_rectangular(tensor);
end

