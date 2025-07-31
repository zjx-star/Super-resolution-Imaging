function mat = rot(x_nodes,y_nodes,z_nodes,ang)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N=x_nodes*y_nodes*z_nodes;
mat=zeros(N,N);
x_mid=ceil(x_nodes/2);
y_mid=ceil(y_nodes/2);
for x=1:x_nodes
    for y=1:y_nodes
        for z=1:z_nodes
            x_loc=x-x_mid;
            y_loc=y-y_mid;
            z_loc=x_loc+y_loc*1i;
            ang_z=angle(z_loc)+ang;
            pho=norm(z_loc);
            x_loc=round(pho*cos(ang_z));
            y_loc=round(pho*sin(ang_z));
            x_rot=x_loc+x_mid;
            y_rot=y_loc+y_mid;
            if(x_rot>=1 && x_rot<=x_nodes && y_rot>=1 && y_rot<=y_nodes)
                n1=(z-1)*x_nodes*y_nodes+(y-1)*x_nodes+x;
                n2=(z-1)*x_nodes*y_nodes+(y_rot-1)*x_nodes+x_rot;
                mat(n2,n1)=1;
            end
        end
    end    
end

end