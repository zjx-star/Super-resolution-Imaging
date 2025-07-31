function accu = asymp_beta_D(m,n,R,h)
% The function is used to estimate beta_mn^D(h). The algorithm 
% is shown in the appendix
%
% Input:
% Precondition:
% Output:
% Postcondition:

m=abs(m);
accu = 1/R*(n*pi/h + (4*m^2-1)*h/(8*n*pi));

end

