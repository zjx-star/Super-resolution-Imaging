function [r,theta,w] = annular_grids(r1,r2,N1,N2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[r,w_g]=lgwt(N1,r1,r2);
theta=[0:N2-1]*2*pi/N2;
w_r=theta*0+2*pi/N2;

[r,theta]=ndgrid(r,theta);
[w_g,w_r]=ndgrid(w_g,w_r);
r=reshape(r,[N1*N2,1]);
theta=reshape(theta,[N1*N2,1]);
w_g=reshape(w_g,[N1*N2,1]);
w_r=reshape(w_r,[N1*N2,1]);
w=w_g.*w_r;
end