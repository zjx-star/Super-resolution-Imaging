function val = H_TM_2(m,n,R,h,l,k,r,theta,mode)
%H_TE_1 computes the first component of H^TE_m0 on the aperture
% 
m1=m;
m=abs(m);
if(mod(mode,2)==0)
    beta=accu_beta_D(m,n,R,h);
    s=s_ij_D(k,m,n,R,h);
    norm=norm_psi_D(m,n,R,h);
    lambda=eigenvalue_D(m,n,R,h);
    val=k/s*2*1i*lambda.^(-1/4).*norm^(-1).*(beta.*(bessely(m,R.*beta).*((besselj(m-1,beta.*r)-besselj(m+1,beta.*r))/2)-...
        besselj(m,R.*beta)*(bessely(m-1,beta.*r)-bessely(m+1,beta.*r))/2).*exp(1i.*m1.*theta).*cos(theta)+...
        1i*m1*(bessely(m,R.*beta).*besselj(m,beta.*r)-besselj(m,R.*beta).*bessely(m,beta.*r)).*exp(1i.*m1.*theta).*(-sin(theta)./r));
else
    beta=accu_beta_D(m,n,R,h);
    s=s_ij_D(k,m,n,R,h);
    norm=norm_psi_D(m,n,R,h);
    lambda=eigenvalue_D(m,n,R,h);
    val=k/s*2*lambda.^(-1/4).*cot(s*l/2)*norm^(-1).*(beta.*(bessely(m,R.*beta).*((besselj(m-1,beta.*r)-besselj(m+1,beta.*r))/2)-...
        besselj(m,R.*beta)*(bessely(m-1,beta.*r)-bessely(m+1,beta.*r))/2).*exp(1i.*m1.*theta).*cos(theta)+...
        1i*m1*(bessely(m,R.*beta).*besselj(m,beta.*r)-besselj(m,R.*beta).*bessely(m,beta.*r)).*exp(1i.*m1.*theta).*(-sin(theta)./r));
end

end

