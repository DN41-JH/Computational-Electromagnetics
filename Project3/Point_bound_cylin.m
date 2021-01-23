% ECE540 PJ3 --- Method of Moment, Jianghuai Liu
% This is the function that returns the spatial information of every
% discretization point along the square boundary contour

function [Xmid_bound_cylin,Ymid_bound_cylin,ds_bound_cylin,angle_bound_cylin]=Point_bound_cylin(Xmin_bound_cylin,Ymin_bound_cylin,Xmax_bound_cylin,Ymax_bound_cylin,Nx_bound_cylin,Ny_bound_cylin)

N_bound_cylin = 2*((Nx_bound_cylin-1)+(Ny_bound_cylin-1));
Xmid_bound_cylin = zeros(N_bound_cylin,1);
Ymid_bound_cylin = zeros(N_bound_cylin,1);
ds_bound_cylin = zeros(N_bound_cylin,1);
angle_bound_cylin = zeros(N_bound_cylin,1);

n = 1;
for I_side = 1:4   
    if (I_side==1)
        Y_temp = linspace(Ymin_bound_cylin,Ymax_bound_cylin,Ny_bound_cylin);
        for j = 1:Ny_bound_cylin-1
            Xmid_bound_cylin(n) = Xmax_bound_cylin;
            Ymid_bound_cylin(n) = (Y_temp(j)+Y_temp(j+1))/2;
            ds_bound_cylin(n) = abs(Y_temp(j)-Y_temp(j+1));
            
            angle_bound_cylin(n) = atan(Ymid_bound_cylin(n)/Xmid_bound_cylin(n));
            angle_bound_cylin(n) = (angle_bound_cylin(n)/pi)*180;
            if (angle_bound_cylin(n)<0)
                angle_bound_cylin(n) = angle_bound_cylin(n)+360;
            end
            n = n+1;
        end
    
    elseif (I_side==2)
        X_temp = linspace(Xmax_bound_cylin,Xmin_bound_cylin,Nx_bound_cylin);
        for i = 1:Nx_bound_cylin-1
            Xmid_bound_cylin(n) = (X_temp(i)+X_temp(i+1))/2;
            Ymid_bound_cylin(n) = Ymax_bound_cylin;
            ds_bound_cylin(n) = abs(X_temp(i)-X_temp(i+1));
            
            angle_bound_cylin(n) = atan(Ymid_bound_cylin(n)/Xmid_bound_cylin(n));
            angle_bound_cylin(n) = (angle_bound_cylin(n)/pi)*180;
            if (angle_bound_cylin(n)<0)
                angle_bound_cylin(n) = angle_bound_cylin(n)+180;
            end
            n = n+1;
        end
        
    elseif (I_side==3)
        Y_temp = linspace(Ymax_bound_cylin,Ymin_bound_cylin,Ny_bound_cylin);
        for j = 1:Ny_bound_cylin-1
            Xmid_bound_cylin(n) = Xmin_bound_cylin;
            Ymid_bound_cylin(n) = (Y_temp(j)+Y_temp(j+1))/2;
            ds_bound_cylin(n) = abs(Y_temp(j)-Y_temp(j+1));
            
            angle_bound_cylin(n) = atan(Ymid_bound_cylin(n)/Xmid_bound_cylin(n));
            angle_bound_cylin(n) = (angle_bound_cylin(n)/pi)*180 + 180;
            n = n+1;
        end        
               
    else
        X_temp = linspace(Xmin_bound_cylin,Xmax_bound_cylin,Nx_bound_cylin);
        for i = 1:Nx_bound_cylin-1
            Xmid_bound_cylin(n) = (X_temp(i)+X_temp(i+1))/2;
            Ymid_bound_cylin(n) = Ymin_bound_cylin;
            ds_bound_cylin(n) = abs(X_temp(i)-X_temp(i+1));
            
            angle_bound_cylin(n) = atan(Ymid_bound_cylin(n)/Xmid_bound_cylin(n));
            angle_bound_cylin(n) = (angle_bound_cylin(n)/pi)*180;
            if (angle_bound_cylin(n)>0)
                angle_bound_cylin(n) = angle_bound_cylin(n)+180;
            else
                angle_bound_cylin(n) = angle_bound_cylin(n)+360;
            end
            n = n+1;
        end
    end
end


