function mat = E3FarField(k,l,X_far,Y_far,Z_far,X_near,Y_near,Z_near,w,E1_near,E2_near,E3_near)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here

N=length(X_far);
M=length(X_near);
mat=zeros(N,M);
for i=1:N
    tar=[X_far(i),Y_far(i),Z_far(i)];
    vec_g13=g13(k,l,tar,X_near,Y_near,Z_near);
    vec_g23=g23(k,l,tar,X_near,Y_near,Z_near);
    vec_g33=g33(k,l,tar,X_near,Y_near,Z_near);
    mat(i,:)=k^2*(E1_near.*vec_g13+E2_near.*vec_g23+E3_near.*vec_g33).*w;
end

end