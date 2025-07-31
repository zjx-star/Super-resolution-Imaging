function mat = trans(x_nodes,y_nodes,z_nodes,trans_x,trans_y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N=x_nodes*y_nodes*z_nodes;
mat=zeros(N,N);
for x=1:x_nodes
    for y=1:y_nodes
        for z=1:z_nodes
            x_trans=x+trans_x;
            y_trans=y+trans_y;
            if(x_trans>=1 && x_trans<=x_nodes && y_trans>=1 && y_trans<=y_nodes)
                n1=(z-1)*x_nodes*y_nodes+(y-1)*x_nodes+x;
                n2=(z-1)*x_nodes*y_nodes+(y_trans-1)*x_nodes+x_trans;
                mat(n2,n1)=1;
            end
        end
    end    
end

end