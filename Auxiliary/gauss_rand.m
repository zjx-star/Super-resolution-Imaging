function vec = gauss_rand(N,std,ratio)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
vec=normrnd(0,std,[N,1])*ratio;
end