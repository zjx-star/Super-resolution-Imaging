function mat = LinearSystem(m,N,l,k,k1,mu_ratio,R1,h1,mode)
%LinearSystem is used to construct the linear system
%Input:
%mode: 0 implies even mode, 1 implies odd mode
%Precondition:
%Output:
%Postcondition
if(m==0)
    %unusable now!!
    mat=zeros(N+1,N+1);
    %mat=zeros(2*N+1,2*N+1);
    mat(1,1)=A_mm(0,l,k,h,mode);
    %mat(1,2:N+1)=R_m_TE(0,N,l,k,h);
    mat(1,2:N+1)=R_m_TM(0,N,1,k,h,mode);
    %mat(2:N+1,1)=C_m_TE(0,N,l,k,h).';
    %mat(2:N+1,2:N+1)=B_m_TE_TE(0,N,l,k,h);
    %mat(2:N+1,N+2:2*N+1)=B_m_TE_TM(0,N,l,k,h);
    mat(2:N+1,1)=C_m_TM(0,N,l,k,h,mode).';
    %mat(N+2:2*N+1,2:N+1)=B_m_TM_TE(0,N,l,k,h);
    mat(2:N+1,2:N+1)=B_m_TM_TM(0,N,l,k,h,mode);
    mat(1:N+2:(N+1)*(N+1))=mat(1:N+2:(N+1)*(N+1))-[D_m(0,l,k,h,mode),ones(1,N)];
else
    if(N~=0)
        mat=zeros((2*N+1),(2*N+1));
        
        mat(1,1)=A_mm(m,l,k,k1,R1,R1,h1,h1,mode);
    
        mat(1,2:N+1)=R_m_TE(m,N,l,k,k1,R1,R1,h1,h1,mode);
        
        mat(1,N+2:2*N+1)=R_m_TM(m,N,l,k,k1,R1,R1,h1,h1,mode);        
                
        mat(2:N+1,1)=C_m_TE(m,N,l,k,k1,R1,R1,h1,h1,mode).';
        
        mat(N+2:2*N+1,1)=C_m_TM(m,N,l,k,k1,R1,R1,h1,h1,mode).';
        
        mat(2:N+1,2:N+1)=B_m_TE_TE(m,N,l,k,k1,R1,R1,h1,h1,mode);
       
        mat(2:N+1,N+2:2*N+1)=B_m_TE_TM(m,N,l,k,k1,R1,R1,h1,h1,mode);
        
        mat(N+2:2*N+1,2:N+1)=B_m_TM_TE(m,N,l,k,k1,R1,R1,h1,h1,mode);
      
        mat(N+2:2*N+1,N+2:2*N+1)=B_m_TM_TM(m,N,l,k,k1,R1,R1,h1,h1,mode);
           
        mat(1:(2*N+2):(2*N+1)*(2*N+1))=mat(1:(2*N+2):(2*N+1)*(2*N+1))-mu_ratio.*[D_m(m,l,k1,R1,h1,mode),ones(1,2*N)];
        mat=-mat;
    else
        %%unusable now
        mat=zeros(1,1);
        mat(1,1)=A_mm(m,l,k,k1,R1,R1,h1,h1,mode);
        mat(1,1)=D_m(m,l,k,k1,R1,h1,mode)-mat(1,1);
    end
end
end