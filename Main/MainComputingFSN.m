I=10;
N=10;
R=0.6*[1:I];
h=0.1./R;
l=1;
mode=0;
d1=sqrt(3)/2; 
k1=[(real(Outm5(1))+real(Outm6(1)))/2,(real(Outm6(1))+real(Outm7(1)))/2,(real(Outm7(1))+real(Outm8(1)))/2,(real(Outm8(1))+real(Outm9(1)))/2,(real(Outm9(1))+real(Outm10(1)))/2];
nk=length(k1);
kinc=k1*0.6;
fourier_series=zeros(I*(2*N+1),nk);
for i=1:nk
    for j=1:I
        RHS=TE_RHS(j,N,l,kinc(i),R(j),h(j),mode,d1);
        fourier_series((j-1)*(2*N+1)+(1:2*N+1),i)=L(:,:,j,i)\RHS;
    end
end