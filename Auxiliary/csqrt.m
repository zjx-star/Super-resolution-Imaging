function val = csqrt(x)
%Compute the square roots of a complex-valued array with the branch defined as (-pi/2,3pi/2]
% 
%Input:
%x: a complex-valued array
%
%Output:
%val: the square roots of x 
%
val=sqrt(x);
Logic=(real(x)<=0)&(imag(x)<0);
val(Logic)=-val(Logic);
end

