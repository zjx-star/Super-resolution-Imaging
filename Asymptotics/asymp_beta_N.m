function val = asymp_beta_N(m,n,R,h)
% The function is used to estimate beta_mn^N(h). The algorithm 
% is shown in the appendix
%
% Input:
% Precondition:
% Output:
% Postcondition:

m=abs(m);

if(n~=0)
    val = 1/R*(n*pi/h + (4*m^2+3)*h/(8*n*pi*(1+h)));
else 
    val = 1/R*(m - m*h/2);
end


end
