function val = H_TE_3(m,n,R,h,l,k,r,theta,mode)
%H_TE_1 computes the first component of H^TE_m0 on the aperture
% 
m1=m;
m=abs(m);
if(mod(mode,2)==0)
    if(n==0)
        beta=accu_beta_N(m,n,R,h);
        s=s_mn_N(k,m,n,R,h);
        lam=eigenvalue_N(m,n,R,h);
        norm=norm_psi_N(m,n,R,h);
        val=1/k*(-1i)*lam*2*cos(s*l/2)*norm^(-1)*(((bessely(m-1,R.*beta)-bessely(m+1,R.*beta))/2.*((besselj(m,beta.*r)))-...
        (besselj(m-1,R.*beta)-besselj(m+1,R.*beta))/2*(bessely(m,beta.*r))).*exp(1i.*m1.*theta));
    else
        beta=accu_beta_N(m,n,R,h);
        s=s_mn_N(k,m,n,R,h);
        lam=eigenvalue_N(m,n,R,h);
        norm=norm_psi_N(m,n,R,h);
        val=1/k*(-1i)*lam^(1/4)*2*cot(s*l/2)*norm^(-1)*(((bessely(m-1,R.*beta)-bessely(m+1,R.*beta))/2.*((besselj(m,beta.*r)))-...
        (besselj(m-1,R.*beta)-besselj(m+1,R.*beta))/2*(bessely(m,beta.*r))).*exp(1i.*m1.*theta));
    end
else
    if(n==0)
        beta=accu_beta_N(m,n,R,h);
        s=s_mn_N(k,m,n,R,h);
        lam=eigenvalue_N(m,n,R,h);
        norm=norm_psi_N(m,n,R,h);
        val=1/k*(-1i)*lam*2i*sin(s*l/2)*norm^(-1)*(((bessely(m-1,R.*beta)-bessely(m+1,R.*beta))/2.*((besselj(m,beta.*r)))-...
        (besselj(m-1,R.*beta)-besselj(m+1,R.*beta))/2*(bessely(m,beta.*r))).*exp(1i.*m1.*theta));
    else
        beta=accu_beta_N(m,n,R,h);
        s=s_mn_N(k,m,n,R,h);
        lam=eigenvalue_N(m,n,R,h);
        norm=norm_psi_N(m,n,R,h);
        val=1/k*(-1i)*lam^(1/4)*2i*norm^(-1)*(((bessely(m-1,R.*beta)-bessely(m+1,R.*beta))/2.*((besselj(m,beta.*r)))-...
        (besselj(m-1,R.*beta)-besselj(m+1,R.*beta))/2*(bessely(m,beta.*r))).*exp(1i.*m1.*theta));        
    end
end

end

