function val = H_TE_1(m,n,R,h,l,k,r,theta,mode)
%H_TE_1 computes the first component of H^TE_m0 on the aperture
% 
m1=m;
m=abs(m);
if(mod(mode,2)==0)
    if(n==0)
        beta=accu_beta_N(m,n,R,h);
        s=s_mn_N(k,m,n,R,h);
        norm=norm_psi_N(m,n,R,h);
        val=s/k*2*1i*sin(s*l/2)*norm^(-1)*(beta.*((bessely(m-1,R.*beta)-bessely(m+1,R.*beta))/2.*((besselj(m-1,beta.*r)-besselj(m+1,beta.*r))/2)-...
        (besselj(m-1,R.*beta)-besselj(m+1,R.*beta))/2*(bessely(m-1,beta.*r)-bessely(m+1,beta.*r))/2).*exp(1i.*m1.*theta).*cos(theta)+...
        1i*m1*((bessely(m-1,R.*beta)-bessely(m+1,R.*beta))/2.*besselj(m,beta.*r)-...
        (besselj(m-1,R.*beta)-besselj(m+1,R.*beta))/2*bessely(m,beta.*r)).*exp(1i.*m1.*theta).*(-sin(theta)./r));
    else
        beta=accu_beta_N(m,n,R,h);
        s=s_mn_N(k,m,n,R,h);
        norm=norm_psi_N(m,n,R,h);
        lambda=eigenvalue_N(m,n,R,h);
        val=s/k*2*1i*lambda.^(-3/4)*norm^(-1)*(beta.*((bessely(m-1,R.*beta)-bessely(m+1,R.*beta))/2.*((besselj(m-1,beta.*r)-besselj(m+1,beta.*r))/2)-...
        (besselj(m-1,R.*beta)-besselj(m+1,R.*beta))/2*(bessely(m-1,beta.*r)-bessely(m+1,beta.*r))/2).*exp(1i.*m1.*theta).*cos(theta)+...
        1i*m1*((bessely(m-1,R.*beta)-bessely(m+1,R.*beta))/2.*besselj(m,beta.*r)-...
        (besselj(m-1,R.*beta)-besselj(m+1,R.*beta))/2*bessely(m,beta.*r)).*exp(1i.*m1.*theta).*(-sin(theta)./r));
    end
else
    if(n==0)
        beta=accu_beta_N(m,n,R,h);
        s=s_mn_N(k,m,n,R,h);
        norm=norm_psi_N(m,n,R,h);
        val=s/k*2*cos(s*l/2)*norm^(-1)*(beta.*((bessely(m-1,R.*beta)-bessely(m+1,R.*beta))/2.*((besselj(m-1,beta.*r)-besselj(m+1,beta.*r))/2)-...
        (besselj(m-1,R.*beta)-besselj(m+1,R.*beta))/2*(bessely(m-1,beta.*r)-bessely(m+1,beta.*r))/2).*exp(1i.*m1.*theta).*cos(theta)+...
        1i*m1*((bessely(m-1,R.*beta)-bessely(m+1,R.*beta))/2.*besselj(m,beta.*r)-...
        (besselj(m-1,R.*beta)-besselj(m+1,R.*beta))/2*bessely(m,beta.*r)).*exp(1i.*m1.*theta).*(-sin(theta)./r));
    else
        beta=accu_beta_N(m,n,R,h);
        s=s_mn_N(k,m,n,R,h);
        norm=norm_psi_N(m,n,R,h);
        lambda=eigenvalue_N(m,n,R,h);
        val=s/k*2*lambda.^(-3/4).*cot(s*l/2)*norm^(-1)*(beta.*((bessely(m-1,R.*beta)-bessely(m+1,R.*beta))/2.*((besselj(m-1,beta.*r)-besselj(m+1,beta.*r))/2)-...
        (besselj(m-1,R.*beta)-besselj(m+1,R.*beta))/2*(bessely(m-1,beta.*r)-bessely(m+1,beta.*r))/2).*exp(1i.*m1.*theta).*cos(theta)+...
        1i*m1*((bessely(m-1,R.*beta)-bessely(m+1,R.*beta))/2.*besselj(m,beta.*r)-...
        (besselj(m-1,R.*beta)-besselj(m+1,R.*beta))/2*bessely(m,beta.*r)).*exp(1i.*m1.*theta).*(-sin(theta)./r));
    end
end

end

