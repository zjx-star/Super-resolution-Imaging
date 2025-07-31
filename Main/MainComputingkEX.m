m=1;                                            %The index of the resonance
n=0;
N=5;                                            %The order of submatrices
R=0.6;
l=1;                                            %Thickness of the slab
h=0.1/R;                                          %The starting h
mode=0;                                         %0 imples even mode while 1 imples odd mode
mu_ratio=1;
eps_ratio=25/9;
delta=1e-6;epsilon=1e-6;max1=10;

k0=kasymptotic(m,n,eps_ratio,mu_ratio,l,R,h,mode);
k1=real(k0);k2=k1+h;k3=k1-h; %Initial values
f2=@(x)(min(eig(LinearSystem(m,N,l,x/sqrt(mu_ratio)/sqrt(eps_ratio),x,mu_ratio,R,h,mode))));
km1=muller(f2,k1,k2,k3,delta,epsilon,max1);
Outm1=[km1,k0,km1-k0];
disp(Outm1)


m=2;
R=0.6*m;
h=0.1/R;

k0=kasymptotic(m,n,eps_ratio,mu_ratio,l,R,h,mode);
k1=real(k0);k2=k1+h;k3=k1-h; %Initial values
f2=@(x)(min(eig(LinearSystem(m,N,l,x/sqrt(mu_ratio)/sqrt(eps_ratio),x,mu_ratio,R,h,mode))));
km2=muller(f2,k1,k2,k3,delta,epsilon,max1);
Outm2=[km2,k0,km2-k0];
disp(Outm2)


m=3;
R=0.6*m;
h=0.1/R;

k0=kasymptotic(m,n,eps_ratio,mu_ratio,l,R,h,mode);
k1=real(k0);k2=k1+h;k3=k1-h; %Initial values
f2=@(x)(min(eig(LinearSystem(m,N,l,x/sqrt(mu_ratio)/sqrt(eps_ratio),x,mu_ratio,R,h,mode))));
km3=muller(f2,k1,k2,k3,delta,epsilon,max1);
Outm3=[km3,k0,km3-k0];
disp(Outm3)


m=4;
R=0.6*m;
h=0.1/R;

k0=kasymptotic(m,n,eps_ratio,mu_ratio,l,R,h,mode);
k1=real(k0);k2=k1+h;k3=k1-h; %Initial values
f2=@(x)(min(eig(LinearSystem(m,N,l,x/sqrt(mu_ratio)/sqrt(eps_ratio),x,mu_ratio,R,h,mode))));
km4=muller(f2,k1,k2,k3,delta,epsilon,max1);
Outm4=[km4,k0,km4-k0];
disp(Outm4)


m=5;
R=0.6*m;
h=0.1/R;

k0=kasymptotic(m,n,eps_ratio,mu_ratio,l,R,h,mode);
k1=real(k0);k2=k1+h;k3=k1-h; %Initial values
f2=@(x)(min(eig(LinearSystem(m,N,l,x/sqrt(mu_ratio)/sqrt(eps_ratio),x,mu_ratio,R,h,mode))));
km5=muller(f2,k1,k2,k3,delta,epsilon,max1);
Outm5=[km5,k0,km5-k0];
disp(Outm5)


m=6;
R=0.6*m;
h=0.1/R;

k0=kasymptotic(m,n,eps_ratio,mu_ratio,l,R,h,mode);
k1=real(k0);k2=k1+h;k3=k1-h; %Initial values
f2=@(x)(min(eig(LinearSystem(m,N,l,x/sqrt(mu_ratio)/sqrt(eps_ratio),x,mu_ratio,R,h,mode))));
km6=muller(f2,k1,k2,k3,delta,epsilon,max1);
Outm6=[km6,k0,km6-k0];
disp(Outm6)


m=7;
R=0.6*m;
h=0.1/R;

k0=kasymptotic(m,n,eps_ratio,mu_ratio,l,R,h,mode);
k1=real(k0);k2=k1+h;k3=k1-h; %Initial values
f2=@(x)(min(eig(LinearSystem(m,N,l,x/sqrt(mu_ratio)/sqrt(eps_ratio),x,mu_ratio,R,h,mode))));
km7=muller(f2,k1,k2,k3,delta,epsilon,max1);
Outm7=[km7,k0,km7-k0];
disp(Outm7)


m=8;
R=0.6*m;
h=0.1/R;

k0=kasymptotic(m,n,eps_ratio,mu_ratio,l,R,h,mode);
k1=real(k0);k2=k1+h;k3=k1-h; %Initial values
f2=@(x)(min(eig(LinearSystem(m,N,l,x/sqrt(mu_ratio)/sqrt(eps_ratio),x,mu_ratio,R,h,mode))));
km8=muller(f2,k1,k2,k3,delta,epsilon,max1);
Outm8=[km8,k0,km8-k0];
disp(Outm8)


m=9;
R=0.6*m;
h=0.1/R;

k0=kasymptotic(m,n,eps_ratio,mu_ratio,l,R,h,mode);
k1=real(k0);k2=k1+h;k3=k1-h; %Initial values
f2=@(x)(min(eig(LinearSystem(m,N,l,x/sqrt(mu_ratio)/sqrt(eps_ratio),x,mu_ratio,R,h,mode))));
km9=muller(f2,k1,k2,k3,delta,epsilon,max1);
Outm9=[km9,k0,km9-k0];
disp(Outm9)


m=10;
R=0.6*m;
h=0.1/R;

k0=kasymptotic(m,n,eps_ratio,mu_ratio,l,R,h,mode);
k1=real(k0);k2=k1+h;k3=k1-h; %Initial values
f2=@(x)(min(eig(LinearSystem(m,N,l,x/sqrt(mu_ratio)/sqrt(eps_ratio),x,mu_ratio,R,h,mode))));
km10=muller(f2,k1,k2,k3,delta,epsilon,max1);
Outm10=[km10,k0,km10-k0];
disp(Outm10)
