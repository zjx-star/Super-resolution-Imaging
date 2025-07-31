function val = norm_psi_D(m,n,R,h)
%norm_psi_D Compute the norm of psi_D eigenfunction
%

beta=accu_beta_D(m,n,R,h);
f=@(x)(bessely(m,beta*R).*besselj(m,beta.*R.*(1+x.*h))-besselj(m,beta*R).*bessely(m,beta.*R.*(1+x.*h))).^2.*R.*(1+x.*h);
I=quad_legendre_seg(f,0,1,n+1);
val=sqrt(2*pi*R*h*I);
end

