function [E1,E2,E3] = plane_wave(E,d,k,X,Y,Z)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E1=E(1)*exp(1i*k*(d(1)*X+d(2)*Y+d(3)*Z));
E2=E(2)*exp(1i*k*(d(1)*X+d(2)*Y+d(3)*Z));
E3=E(3)*exp(1i*k*(d(1)*X+d(2)*Y+d(3)*Z));

end