function mat = LinearSystemN(m,N,I,l,k,k1,mu_ratio,R,h,mode)
%LinearSystem is used to construct the linear system
%Input:
%mode: 0 implies even mode, 1 implies odd mode
%Precondition:
%Output:
%Postcondition
Nu=2*N+1;
mat=zeros(I*Nu,I*Nu);

for i=1:I
    for j=1:I
        mat((i-1)*Nu+1,(j-1)*Nu+1)=A_mm(m,l,k,k1,R(i),R(j),h(i),h(j),mode);

        mat((i-1)*Nu+1,(j-1)*Nu+(2:N+1))=R_m_TE(m,N,l,k,k1,R(i),R(j),h(i),h(j),mode);
        
        mat((i-1)*Nu+1,(i-1)*Nu+(N+2:2*N+1))=R_m_TM(m,N,l,k,k1,R(i),R(j),h(i),h(j),mode);        
                
        mat((i-1)*Nu+(2:N+1),1)=C_m_TE(m,N,l,k,k1,R(i),R(j),h(i),h(j),mode).';
        
        mat((i-1)*Nu+(N+2:2*N+1),1)=C_m_TM(m,N,l,k,k1,R(i),R(j),h(i),h(j),mode).';
        
        mat((i-1)*Nu+(2:N+1),(j-1)*Nu+(2:N+1))=B_m_TE_TE(m,N,l,k,k1,R(i),R(j),h(i),h(j),mode);
       
        mat((i-1)*Nu+(2:N+1),(j-1)*Nu+(N+2:2*N+1))=B_m_TE_TM(m,N,l,k,k1,R(i),R(j),h(i),h(j),mode);
        
        mat((i-1)*Nu+(N+2:2*N+1),(j-1)*Nu+(2:N+1))=B_m_TM_TE(m,N,l,k,k1,R(i),R(j),h(i),h(j),mode);
      
        mat((i-1)*Nu+(N+2:2*N+1),(j-1)*Nu+(N+2:2*N+1))=B_m_TM_TM(m,N,l,k,k1,R(i),R(j),h(i),h(j),mode);
        if(i==j)
            mat((i-1)*Nu+(1:(2*N+1)),(j-1)*Nu+(1:(2*N+1)))=mat((i-1)*Nu+(1:(2*N+1)),(j-1)*Nu+(1:(2*N+1)))-mu_ratio.*diag([D_m(m,l,k1,R(i),h(i),mode),ones(1,2*N)]);
        end
    end
end
        mat=-mat;
end
