function mat = s_mn_N(k,m,n,R,h)
%Compute s_mn^N defined in Theorem 2.4
%
%Input: 
%k: wavenumber
%m: angular momentum of the eigenvalue defined in eq. (3.11)
%n: the longitudinal index of the eigenvalue defined in eq. (3.11)
%R,h: the eigenvalue is in an annular aperture of radius R and width R*h
%
%Precondition:
%
%Output:
%mat: the value of s_mn^N
%
%Postcondition:
%
m=abs(m);
eig_mat=eigenvalue_N(m,n,R,h);
mat=csqrt(k^2-eig_mat);

end

