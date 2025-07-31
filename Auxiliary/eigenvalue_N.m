function mat = eigenvalue_N(m,n,R,h,delta,epsilon,maxi)
%Compute the eigenvalue defined in eq. (3.11)
%
%Input: 
%m: angular momentum of the eigenvalue defined in eq. (3.11)
%n: the longitudinal index of the eigenvalue defined in eq. (3.11)
%R,h: the eigenvalue is in an annular aperture of radius R and width R*h
%delta,epsilon,maxi: stop criteria for computing beta_N by Newton's method
%
%Precondition: 
%m,n are all row vectors
%
%Output:
%mat: the eigenvalue defined in eq. (3.11)
%
%Postcondition:
%

if(nargin<5)
    delta=1e-10;
    epsilon=1e-10;
    maxi=100;
end

m=abs(m);

l1=length(m);
l2=length(n);
[X,Y]=meshgrid(m,n);
mat=zeros(l2,l1);
for j=[1:1:l2]
    for i=[1:1:l1]
        mat(j,i)=(accu_beta_N(abs(X(j,i)),Y(j,i),R,h,delta,epsilon,maxi))^2;
    end
end
mat=mat';
end

