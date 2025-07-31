%% imaging settings
l=1;
I=10;
Nmodes=10;
R=0.6*[1:I];
k=[(real(Outm5(1))+real(Outm6(1)))/2,(real(Outm6(1))+real(Outm7(1)))/2,(real(Outm7(1))+real(Outm8(1)))/2,(real(Outm8(1))+real(Outm9(1)))/2,(real(Outm9(1))+real(Outm10(1)))/2];
nk=length(k);
kinc=0.6*k;
h=0.1./R;
detector_height=20*pi;
Nodes=120;
N_far=(Nodes+1)^2;
detector_range=[-12*pi:24*pi/Nodes:12*pi];
[X_far,Y_far]=ndgrid(detector_range,detector_range);
Z_far=0*X_far+l/2+detector_height;
X_far=reshape(X_far,[1,N_far]);
Y_far=reshape(Y_far,[1,N_far]);
Z_far=reshape(Z_far,[1,N_far]);

[X_near,Y_near,Z_near,w_near]=Grid2([-sqrt(2)*pi,sqrt(2)*pi-pi,l/2+0.1],[-sqrt(2)*pi+pi,sqrt(2)*pi,l/2+0.3],40,40,4);
N_near=40*40*4;

%% near field computing
% E1_near_k=zeros(nk,N_near);
% E2_near_k=zeros(nk,N_near);
% E3_near_k=zeros(nk,N_near);
% 
% for i=1:nk
%     temp=E1NearFieldN_Grid(kinc(i),k(i),R,h,l,X_near,Y_near,Z_near,10,I,fourier_series(:,i));
%     E1_near_k(i,:)=temp;
%     temp=E2NearFieldN_Grid(kinc(i),k(i),R,h,l,X_near,Y_near,Z_near,10,I,fourier_series(:,i));
%     E2_near_k(i,:)=temp;
%     temp=E3NearFieldN_Grid(kinc(i),k(i),R,h,l,X_near,Y_near,Z_near,10,I,fourier_series(:,i));
%     E3_near_k(i,:)=temp;
% end

%%plane waves
% for i=1:nk
% [E1,E2,E3]=plane_wave([1,0,0],[0,1,0],kinc(i),X_near,Y_near,Z_near);
% E1_near_k(i,:)=E1;
% E2_near_k(i,:)=E2;
% E3_near_k(i,:)=E3;
% end


%% far field matrix caused by source
% E1_far_mat_k=zeros(N_far,N_near,nk);
% E2_far_mat_k=zeros(N_far,N_near,nk);
% E3_far_mat_k=zeros(N_far,N_near,nk);
% 
% for i=1:nk
%     E1_far_mat_k(:,:,i)=E1FarField(kinc(i),l,X_far,Y_far,Z_far,X_near,Y_near,Z_near,w_near,E1_near_k(i,:),E2_near_k(i,:),E3_near_k(i,:));
%     E2_far_mat_k(:,:,i)=E2FarField(kinc(i),l,X_far,Y_far,Z_far,X_near,Y_near,Z_near,w_near,E1_near_k(i,:),E2_near_k(i,:),E3_near_k(i,:));
%     E3_far_mat_k(:,:,i)=E3FarField(kinc(i),l,X_far,Y_far,Z_far,X_near,Y_near,Z_near,w_near,E1_near_k(i,:),E2_near_k(i,:),E3_near_k(i,:));  
% end


%% Grids on the aperture
% [r_aper,theta_aper,w_aper]=annular_grids(0,0.1,8,48);
% N_aper=8*48;
% R_aper=zeros(N_aper*I,1);
% Theta_aper=zeros(N_aper*I,1);
% W_aper=zeros(N_aper*I,1);
% for i=1:I
%     R_aper((i-1)*N_aper+(1:N_aper))=r_aper+R(i);
%     Theta_aper((i-1)*N_aper+(1:N_aper))=theta_aper;
%     W_aper((i-1)*N_aper+(1:N_aper))=w_aper;
% end
% X_aper=R_aper.*cos(Theta_aper);
% Y_aper=R_aper.*sin(Theta_aper);
% Z_aper=X_aper*0+l/2;


%% H field on the aperture 
% H1_aper_mat_k=zeros(I*N_aper,N_near,nk);
% H2_aper_mat_k=zeros(I*N_aper,N_near,nk);
% 
% for i=1:nk
%     temp=H1FarField(kinc(i),l,X_aper,Y_aper,Z_aper,X_near,Y_near,Z_near,w_near,E1_near_k(i,:),E2_near_k(i,:),E3_near_k(i,:));
%     H1_aper_mat_k(:,:,i)=2*temp;
%     temp=H2FarField(kinc(i),l,X_aper,Y_aper,Z_aper,X_near,Y_near,Z_near,w_near,E1_near_k(i,:),E2_near_k(i,:),E3_near_k(i,:));
%     H2_aper_mat_k(:,:,i)=2*temp;
% end


%% resoant fourier_series matrix
% fourier_coe_mat=zeros(I*(2*Nmodes+1),N_near,nk);
% for j=1:nk
% for i=1:I
%     lambda=eigenvalue_N(i,0,R(i),h(i));
%     [par1_Psi,par2_Psi]=nabla_Psi_N(i,0,R(i),h(i),R_aper((i-1)*N_aper+(1:N_aper)),Theta_aper((i-1)*N_aper+(1:N_aper)));
%     fourier_coe_mat((i-1)*(2*Nmodes+1)+1,:,j)=-1/2/lambda*(((conj(par1_Psi).*R_aper((i-1)*N_aper+(1:N_aper)).*W_aper((i-1)*N_aper+(1:N_aper))).')*H1_aper_mat_k((i-1)*N_aper+(1:N_aper),:,j)+...
%         ((conj(par2_Psi).*R_aper((i-1)*N_aper+(1:N_aper)).*W_aper((i-1)*N_aper+(1:N_aper))).')*H2_aper_mat_k((i-1)*N_aper+(1:N_aper),:,j));
%     fourier_coe_mat((i-1)*(2*Nmodes+1)+(1:2*Nmodes+1),:,j)=L(:,:,i,j)\fourier_coe_mat((i-1)*(2*Nmodes+1)+(1:2*Nmodes+1),:,j);
% end
% end
% 
% E1_far_res_fourier_mat_k=zeros(N_far,I*(2*Nmodes+1),nk);
% E2_far_res_fourier_mat_k=zeros(N_far,I*(2*Nmodes+1),nk);
% E3_far_res_fourier_mat_k=zeros(N_far,I*(2*Nmodes+1),nk);
% 
% for i=1:nk
%     E1_far_res_fourier_mat_k(:,:,i)=E1FarField_fourier_mat(kinc(i),k(i),l,X_far,Y_far,Z_far,R,h);
%     E2_far_res_fourier_mat_k(:,:,i)=E2FarField_fourier_mat(kinc(i),k(i),l,X_far,Y_far,Z_far,R,h);
%     E3_far_res_fourier_mat_k(:,:,i)=E3FarField_fourier_mat(kinc(i),k(i),l,X_far,Y_far,Z_far,R,h);
% end
% 
% E1_far_res_mat_k=zeros(N_far,N_near,nk);
% E2_far_res_mat_k=zeros(N_far,N_near,nk);
% E3_far_res_mat_k=zeros(N_far,N_near,nk);
% 
% for i=1:nk
%     E1_far_res_mat_k(:,:,i)=E1_far_res_fourier_mat_k(:,:,i)*fourier_coe_mat(:,:,i);
%     E2_far_res_mat_k(:,:,i)=E2_far_res_fourier_mat_k(:,:,i)*fourier_coe_mat(:,:,i);
%     E3_far_res_mat_k(:,:,i)=E3_far_res_fourier_mat_k(:,:,i)*fourier_coe_mat(:,:,i);
% end


%% sum matrix
E1_far_sum_mat_k=zeros(N_far,N_near,nk);
E2_far_sum_mat_k=zeros(N_far,N_near,nk);
E3_far_sum_mat_k=zeros(N_far,N_near,nk);

for i=1:nk
    E1_far_sum_mat_k(:,:,i)=E1_far_mat_k(:,:,i);%+E1_far_res_mat_k(:,:,i);
    E2_far_sum_mat_k(:,:,i)=E2_far_mat_k(:,:,i);%+E2_far_res_mat_k(:,:,i);
    E3_far_sum_mat_k(:,:,i)=E3_far_mat_k(:,:,i);%+E3_far_res_mat_k(:,:,i);       
end


%% gradient decent
p_real=pfunc(0.5,X_near,Y_near,Z_near,2*pi/10,0).';
p0=p_real*0;
A=zeros(N_near,N_near);
for i=1:nk
    A=A+real(E1_far_sum_mat_k(:,:,i))'*real(E1_far_sum_mat_k(:,:,i))+imag(E1_far_sum_mat_k(:,:,i))'*imag(E1_far_sum_mat_k(:,:,i))+...
        real(E2_far_sum_mat_k(:,:,i))'*real(E2_far_sum_mat_k(:,:,i))+imag(E2_far_sum_mat_k(:,:,i))'*imag(E2_far_sum_mat_k(:,:,i))+...
        real(E3_far_sum_mat_k(:,:,i))'*real(E3_far_sum_mat_k(:,:,i))+imag(E3_far_sum_mat_k(:,:,i))'*imag(E3_far_sum_mat_k(:,:,i));
end
rot1=rot(40,40,4,pi/6);
rot2=rot(40,40,4,2*pi/6);
rot3=rot(40,40,4,3*pi/6);
rot4=rot(40,40,4,4*pi/6);
rot5=rot(40,40,4,5*pi/6);
rot6=rot(40,40,4,pi);
rot7=rot(40,40,4,7*pi/6);
rot8=rot(40,40,4,8*pi/6);
rot9=rot(40,40,4,9*pi/6);
rot10=rot(40,40,4,10*pi/6);
rot11=rot(40,40,4,11*pi/6);
A_rot=A+rot1'*A*rot1+rot2'*A*rot2+rot3'*A*rot3+rot4'*A*rot4+rot5'*A*rot5+rot6'*A*rot6...
    +rot7'*A*rot7+rot8'*A*rot8+rot9'*A*rot9+rot10'*A*rot10+rot11'*A*rot11;

% tran1=trans(40,40,4,4,-4);
% tran2=trans(40,40,4,4,4);
% tran3=trans(40,40,4,-4,4);
% tran4=trans(40,40,4,-4,-4);
% 
% tran5=trans(40,40,4,8,-8);
% tran6=trans(40,40,4,8,8);
% tran7=trans(40,40,4,-8,8);
% tran8=trans(40,40,4,-8,-8);

%A1=A_rot;%+tran1.'*A_rot*tran1;%+tran3.'*A_rot*tran3;
A1=A_rot;
%+tran2.'*A_rot*tran2+tran3.'*A_rot*tran3+tran4.'*A_rot*tran4;
   %tran5.'*A_rot*tran5+tran6.'*A_rot*tran6+tran7.'*A_rot*tran7+tran8.'*A_rot*tran8;   
g_real=zeros(N_near,1);
ratio=0.1;
n_rot=12;
for j=1:n_rot
for i=1:nk
    rot1=rot(40,40,4,(j-1)*pi/6);
    g1=real(E1_far_sum_mat_k(:,:,i))*rot1*p_real;
    std1=std(g1);
    g_real=g_real+rot1'*real(E1_far_sum_mat_k(:,:,i))'*(g1+gauss_rand(N_far,std1,ratio));

    g2=imag(E1_far_sum_mat_k(:,:,i))*rot1*p_real;
    std2=std(g2);
    g_real=g_real+rot1'*imag(E1_far_sum_mat_k(:,:,i))'*(g2+gauss_rand(N_far,std2,ratio));

    g3=real(E2_far_sum_mat_k(:,:,i))*rot1*p_real;
    std3=std(g3);
    g_real=g_real+rot1'*real(E2_far_sum_mat_k(:,:,i))'*(g3+gauss_rand(N_far,std3,ratio)); 

    g4=imag(E2_far_sum_mat_k(:,:,i))*rot1*p_real;
    std4=std(g4);
    g_real=g_real+rot1'*imag(E2_far_sum_mat_k(:,:,i))'*(g4+gauss_rand(N_far,std4,ratio));

    g5=real(E3_far_sum_mat_k(:,:,i))*rot1*p_real;
    std5=std(g5);
    g_real=g_real+rot1'*real(E3_far_sum_mat_k(:,:,i))'*(g5+gauss_rand(N_far,std5,ratio)); 

    g6=imag(E3_far_sum_mat_k(:,:,i))*rot1*p_real;
    std6=std(g6);
    g_real=g_real+rot1'*imag(E3_far_sum_mat_k(:,:,i))'*(g6+gauss_rand(N_far,std6,ratio));
end
end

% A_rand=zeros(N_near,N_near);
% for i=1:nk
%     rand1=gauss_rand(N_far,1,ratio);
%     rand2=gauss_rand(N_far,1,ratio);
%     rand3=gauss_rand(N_far,1,ratio);
%     rand4=gauss_rand(N_far,1,ratio);
%     rand5=gauss_rand(N_far,1,ratio);
%     rand6=gauss_rand(N_far,1,ratio);
%     A_rand=A_rand+real(E1_far_sum_mat_k(:,:,i))'*diag(rand1+1)*real(E1_far_sum_mat_k(:,:,i))+imag(E1_far_sum_mat_k(:,:,i))'*diag(rand2+1)*imag(E1_far_sum_mat_k(:,:,i))+...
%         real(E2_far_sum_mat_k(:,:,i))'*diag(rand3+1)*real(E2_far_sum_mat_k(:,:,i))+imag(E2_far_sum_mat_k(:,:,i))'*diag(rand4+1)*imag(E2_far_sum_mat_k(:,:,i))+...
%         real(E3_far_sum_mat_k(:,:,i))'*diag(rand5+1)*real(E3_far_sum_mat_k(:,:,i))+imag(E3_far_sum_mat_k(:,:,i))'*diag(rand6+1)*imag(E3_far_sum_mat_k(:,:,i));
% end
% A_rot_rand=A_rand+rot1'*A_rand*rot1+rot2'*A_rand*rot2+rot3'*A_rand*rot3+rot4'*A_rand*rot4+rot5'*A_rand*rot5+rot6'*A_rand*rot6...
%     +rot7'*A_rand*rot7+rot8'*A_rand*rot8+rot9'*A_rand*rot9+rot10'*A_rand*rot10+rot11'*A_rand*rot11;
% g_real=A_rot_rand*p_real;
% 
for(j=1:5000)
    r=g_real-A1*p0;
    alpha=r.'*r/((A1*r).'*r);
    p1=p0+alpha*r;
    rel_err=norm(alpha*r)/norm(p1);
    display(rel_err);
    if(rel_err < 1e-3)
        break
    end
    p0=p1;
end
