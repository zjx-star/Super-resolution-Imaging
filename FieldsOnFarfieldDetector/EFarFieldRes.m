function [E1_mat,E2_mat,E3_mat] = EFarFieldRes(Nholes,Nmat,Lk,kinc,kma,l,X_far,Y_far,Z_far,R,h,R_aper,Theta_aper,W_aper,H1_aper,H2_aper)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
s=size(H1_aper);
N_near=s(2);
N_aper=s(1)/Nholes;
N_far=length(X_far);
fourier_coe=zeros(Nholes*(2*Nmat+1),N_near);
for i=1:Nholes
    lambda=eigenvalue_N(i,0,R(i),h(i));
    [par1_Psi,par2_Psi]=nabla_Psi_N(i,0,R(i),h(i),R_aper((i-1)*N_aper+(1:N_aper)),Theta_aper((i-1)*N_aper+(1:N_aper)));
    fourier_coe((i-1)*(2*Nmat+1)+1,:)=-1/2/lambda*(((conj(par1_Psi).*R_aper((i-1)*N_aper+(1:N_aper)).*W_aper((i-1)*N_aper+(1:N_aper))).')*H1_aper((i-1)*N_aper+(1:N_aper),:)+...
        ((conj(par2_Psi).*R_aper((i-1)*N_aper+(1:N_aper)).*W_aper((i-1)*N_aper+(1:N_aper))).')*H2_aper((i-1)*N_aper+(1:N_aper),:));
    fourier_coe((i-1)*(2*Nmat+1)+(1:2*Nmat+1),:)=Lk(:,:,i)\fourier_coe((i-1)*(2*Nmat+1)+(1:2*Nmat+1),:);
end

E1_mat=zeros(N_far,N_near);
E2_mat=zeros(N_far,N_near);
E3_mat=zeros(N_far,N_near);
for i=1:N_near
    E1_mat(:,i)=E1NearFieldN_Grid(kinc,kma,R,h,l,X_far,Y_far,Z_far,Nmat,Nholes,fourier_coe(:,i));
    E2_mat(:,i)=E2NearFieldN_Grid(kinc,kma,R,h,l,X_far,Y_far,Z_far,Nmat,Nholes,fourier_coe(:,i));
    E3_mat(:,i)=E3NearFieldN_Grid(kinc,kma,R,h,l,X_far,Y_far,Z_far,Nmat,Nholes,fourier_coe(:,i));    
end
    
end