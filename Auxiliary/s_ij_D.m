function mat = s_ij_D(k,i,j,R,h)
%Compute s_ij^D defined in Theorem 2.4
%
%Input: 
%k: wavenumber
%i: angular momentum of the eigenvalue defined in eq. (3.9)
%j: the longitudinal index of the eigenvalue defined in eq. (3.9)
%R,h: the eigenvalue is in an annular aperture of radius R and width R*h
%
%Precondition:
%
%Output:
%mat: the value of s_mn^N
%
%Postcondition:
%
i=abs(i);
eig_mat=eigenvalue_D(i,j,R,h);
mat=csqrt(k^2-eig_mat);

end

