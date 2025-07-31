function mat = D_m(m,l,k,R,h,mode)
%Construct D_m
%   
%Input:
%mode: 0 imples even mode while 1 imples odd mode
%Precondition:
%Output
%Postcondition

if(m==0)
    mat=sin(k*l/2+pi/2*mode);
else
    s=s_mn_N(k,m,0,R,h);
    mat=s*sin(s*l/2+pi/2*mode);
end

end

