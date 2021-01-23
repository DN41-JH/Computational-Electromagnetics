% ECE540 PJ3 --- Method of Moment, Jianghuai Liu
c = 3*10^8;
Ima_Unit = sqrt(-1);
lambda = 0.01; % The wavelength
k = 2*pi/lambda; % The wavenumber
Z0 = 376.73; % Characteristic impedance of free-space
gamma = 1.781; % Constant in Small-argument approximation of Hankel function
Distance_Far_Field = 20*lambda;
N_angle_far_field = 721; % Number of sampling points of angle around far-field cirlce
angle_far_field = linspace(0,2*pi,N_angle_far_field);

Xmin_domain = -5*lambda; % Computational Domain Size
Ymin_domain = -5*lambda; % Computational Domain Size
Xmax_domain = 5*lambda; % Computational Domain Size
Ymax_domain = 5*lambda; % Computational Domain Size
Nx_domain = 201; % Spatial Discretization of Computational Domain
Ny_domain = 201; % Spatial Discretization of Computational Domain
dx_domain = (Xmax_domain-Xmin_domain)/(Nx_domain-1);
dy_domain = (Ymax_domain-Ymin_domain)/(Ny_domain-1);
X_domain = linspace(Xmin_domain,Xmax_domain,Nx_domain);
Y_domain = linspace(Ymax_domain,Ymin_domain,Ny_domain);

Xlength_bound_cylin = 1*lambda; % Length of Conducting Square Cylinder
Ylength_bound_cylin = 1*lambda; % Width of Conducting Square Cylinder
Xmin_bound_cylin = -0.5*lambda;
Ymin_bound_cylin = -0.5*lambda;
Xmax_bound_cylin = Xmin_bound_cylin + Xlength_bound_cylin;
Ymax_bound_cylin = Ymin_bound_cylin + Ylength_bound_cylin;
Nx_bound_cylin = 101; % Spatial Discretization of Boundary Contour of Square Cylinder
Ny_bound_cylin = 101; % Spatial Discretization of Boundary Contour of Square Cylinder
N_bound_cylin = 2*((Nx_bound_cylin-1)+(Ny_bound_cylin-1));

[Xmid_bound_cylin,Ymid_bound_cylin,ds_bound_cylin,angle_bound_cylin]=Point_bound_cylin(Xmin_bound_cylin,Ymin_bound_cylin,Xmax_bound_cylin,Ymax_bound_cylin,Nx_bound_cylin,Ny_bound_cylin);

TE = 1;
TM = 0;

Z = zeros(N_bound_cylin);
V = zeros(N_bound_cylin,1);

if (TM==1)
    Ez_sc = zeros(Nx_domain,Ny_domain);
    Ez_inc = zeros(Nx_domain,Ny_domain);
    Ez_total = zeros(Nx_domain,Ny_domain);
    Ez_sc_far = zeros(N_angle_far_field,1);
    
    for m = 1:N_bound_cylin
        for n = 1:N_bound_cylin
            if (m==n)
                Z(m,n) = (k*Z0*ds_bound_cylin(n)/4)*(1-(2i/pi)*log(k*gamma*ds_bound_cylin(n)/(4*exp(1))));
            else
                distance = sqrt((Xmid_bound_cylin(m)-Xmid_bound_cylin(n))^2 + (Ymid_bound_cylin(m)-Ymid_bound_cylin(n))^2);
                Z(m,n) = (k*Z0*ds_bound_cylin(n)/4)*besselh(0,2,k*distance);
            end
        V(m) = exp(-Ima_Unit*(k*Xmid_bound_cylin(m))); % Incident Plane Wave propagating in -(x-hat) direction
        %V(m) = exp(-Ima_Unit*k*(Xmid_bound_cylin(m)+Ymid_bound_cylin(m))/sqrt(2)); % Incident Plane Wave propagating in -[(x-hat + y-hat)/sqrt(2)] direction
        end
    end
    
    Jz = linsolve(Z,V); % Solve For the Induced Surface Current Density
       
    shift_count = 1; % Countour Angle Shifting Counter
    for n = 1:N_bound_cylin
        if (angle_bound_cylin(n)>min(angle_bound_cylin))
            shift_count = shift_count + 1; 
        else
            shift_count = shift_count - 1;
            break
        end
    end
    
    angle_bound_cylin = circshift(angle_bound_cylin,-1*shift_count);
    Xmid_bound_cylin = circshift(Xmid_bound_cylin,-1*shift_count);
    Ymid_bound_cylin = circshift(Ymid_bound_cylin,-1*shift_count);
    ds_bound_cylin = circshift(ds_bound_cylin,-1*shift_count);
    Jz = circshift(Jz,-1*shift_count);

    
    for I_far = 1:N_angle_far_field % Calculating the Far-Field Scattered Field
        x_far_field = Distance_Far_Field*cos(angle_far_field(I_far));
        y_far_field = Distance_Far_Field*sin(angle_far_field(I_far));
        for s = 1:N_bound_cylin
            distance = sqrt((x_far_field-Xmid_bound_cylin(s))^2 + (y_far_field-Ymid_bound_cylin(s))^2);
            Ez_sc_far(I_far) = Ez_sc_far(I_far) - (k*Z0*ds_bound_cylin(s)/4)*besselh(0,2,k*distance)*Jz(s);
        end
    end

    
    for i = 1:Nx_domain % Calculating the Scattered Field and Total Field on the Entire Computational Domain
        for j = 1:Ny_domain
            if ((X_domain(i)>Xmin_bound_cylin)&&(X_domain(i)<Xmax_bound_cylin)&&(Y_domain(j)>Ymin_bound_cylin)&&(Y_domain(j)<Ymax_bound_cylin))
                Ez_sc(i,j) = 0;
                Ez_inc(i,j) = 0;
            else
                for s = 1:N_bound_cylin
                    distance = sqrt((X_domain(i)-Xmid_bound_cylin(s))^2 + (Y_domain(j)-Ymid_bound_cylin(s))^2);
                    Ez_sc(i,j) = Ez_sc(i,j) - (k*Z0*ds_bound_cylin(s)/4)*besselh(0,2,k*distance)*Jz(s);
                end
                Ez_inc(i,j) = exp(-Ima_Unit*k*X_domain(i)); % Incident Plane Wave propagating in -(x-hat) direction
                %Ez_inc(i,j) = exp(-Ima_Unit*k*(X_domain(i)+Y_domain(j))/sqrt(2)); % Incident Plane Wave propagating in -[(x-hat + y-hat)/sqrt(2)] direction
                Ez_total(i,j) = Ez_inc(i,j) + Ez_sc(i,j);
            end
        end
    end
    
    % Next Goes For Plotting Commands
    
    %plot(angle_bound_cylin,abs(Jz));
    %plot(angle_bound_cylin,abs(Ez_sc_far));
    %xlim([0 360]);
    %xlabel('Azimuthal Angle \phi   [\circ]','Fontsize',14);
    %ylabel('Induced Surface Current Density |J_z|    [A/m^2]','Fontsize',14);
    %title('Under Incident Plane Wave Propagating in $-\hat{x}$ (TM Polarization)  [ECE540 PJ3---Jianghuai Liu]','Interpreter','Latex','Fontsize',18);
    %title('Under Incident Plane Wave Propagating in $-\frac{\hat{x}+\hat{y}}{\sqrt{2}}$ (TM Polarization)  [ECE540 PJ3---Jianghuai Liu]','Interpreter','Latex','Fontsize',18);
    
    
    polarplot(angle_far_field,circshift(abs(Ez_sc_far),(N_angle_far_field-1)/2));
    title('Far Field (at r=20$\lambda$) Scattered $|{E_z^{sc}}|$ Under Incident Plane Wave Propagating in $-\hat{x}$ (TM Polarization)  [ECE540 PJ3---Jianghuai Liu]','Interpreter','Latex','Fontsize',18);
    %title('Far Field (at r=20$\lambda$) Scattered $|{E_z^{sc}}|$ Under Incident Plane Wave Propagating in $-\frac{\hat{x}+\hat{y}}{\sqrt{2}}$ (TM Polarization)  [ECE540 PJ3---Jianghuai Liu]','Interpreter','Latex','Fontsize',18);
    
    
    %[xPlot,yPlot] = meshgrid(X_domain./lambda,Y_domain./lambda);
    %contour(yPlot,xPlot,abs(Ez_total),6000);
    %xlabel('X/\lambda');
    %ylabel('Y/\lambda');
    %title('Total Electric Field $\left|E_z\right|$ in [V/m], with Incident Plane Wave Propagating in $-\hat{x}$ (TM Polarization)  [ECE540 PJ3---Jianghuai Liu]','Interpreter','Latex','Fontsize',14);
    %title('Total Electric Field $\left|E_z\right|$ in [V/m], with Incident Plane Wave Propagating in $-\frac{\hat{x}+\hat{y}}{\sqrt{2}}$ (TM Polarization)  [ECE540 PJ3---Jianghuai Liu]','Interpreter','Latex','Fontsize',14);
    %colorbar();
end


if (TE==1)
    Hz_sc = zeros(Nx_domain,Ny_domain);
    Hz_inc = zeros(Nx_domain,Ny_domain);
    Hz_total = zeros(Nx_domain,Ny_domain);
    Hz_sc_far = zeros(N_angle_far_field,1);
    
    for m = 1:N_bound_cylin
        for n = 1:N_bound_cylin
            if (m==n)
                Z(m,n) = -0.5;
            else
                distance = sqrt((Xmid_bound_cylin(m)-Xmid_bound_cylin(n))^2 + (Ymid_bound_cylin(m)-Ymid_bound_cylin(n))^2);
                if ((n>=1) && (n<=(Ny_bound_cylin-1))) % First segment, normal vector is x-hat
                    Z(m,n) = (k*ds_bound_cylin(n)/(4i))*besselh(1,2,k*distance)*(Xmid_bound_cylin(m)-Xmid_bound_cylin(n))/distance;
                elseif ((n>=Ny_bound_cylin) && (n<=(Nx_bound_cylin+Ny_bound_cylin-2))) % Second segment, normal vector is y-hat
                    Z(m,n) = (k*ds_bound_cylin(n)/(4i))*besselh(1,2,k*distance)*(Ymid_bound_cylin(m)-Ymid_bound_cylin(n))/distance;
                elseif ((n>=(Nx_bound_cylin+Ny_bound_cylin-2)) && (n<=(Nx_bound_cylin+2*Ny_bound_cylin-3))) % Third segment, normal vector is -x-hat
                    Z(m,n) = (k*ds_bound_cylin(n)/(4i))*besselh(1,2,k*distance)*(Xmid_bound_cylin(n)-Xmid_bound_cylin(m))/distance;
                else % Fourth segment, normal vector is -y-hat
                    Z(m,n) = (k*ds_bound_cylin(n)/(4i))*besselh(1,2,k*distance)*(Ymid_bound_cylin(n)-Ymid_bound_cylin(m))/distance;
                end
            end
        %V(m) = exp(-Ima_Unit*k*Xmid_bound_cylin(m)); % Incident Plane Wave Propagating in -(x-hat) direction
        V(m) = exp(-Ima_Unit*k*(Xmid_bound_cylin(m)+Ymid_bound_cylin(m))/sqrt(2)); % Incident Plane Wave propagating in -[(x-hat + y-hat)/sqrt(2)] direction
        end
    end
    
    Jt = linsolve(Z,V); % Solve for the Induced Surface Current Density
    
    for I_far = 1:N_angle_far_field % Calculating the Far-Field Scattered Field On the Far-Field Circle
        x_far_field = Distance_Far_Field*cos(angle_far_field(I_far));
        y_far_field = Distance_Far_Field*sin(angle_far_field(I_far));
        
        for s = 1:N_bound_cylin
            distance = sqrt((x_far_field-Xmid_bound_cylin(s))^2 + (y_far_field-Ymid_bound_cylin(s))^2);
            if ((s>=1) && (s<=(Ny_bound_cylin-1))) % First segment, normal vector is x-hat
                Hz_sc_far(I_far) = Hz_sc_far(I_far) - (k*ds_bound_cylin(s)/(4i))*besselh(1,2,k*distance)*Jt(s)*(x_far_field-Xmid_bound_cylin(s))/distance;
            elseif ((s>=Ny_bound_cylin) && (s<=(Nx_bound_cylin+Ny_bound_cylin-2))) % Second segment, normal vector is y-hat
                Hz_sc_far(I_far) = Hz_sc_far(I_far) - (k*ds_bound_cylin(s)/(4i))*besselh(1,2,k*distance)*Jt(s)*(y_far_field-Ymid_bound_cylin(s))/distance;
            elseif ((s>=(Nx_bound_cylin+Ny_bound_cylin-2)) && (s<=(Nx_bound_cylin+2*Ny_bound_cylin-3))) % Third segment, normal vector is -x-hat
                Hz_sc_far(I_far) = Hz_sc_far(I_far) - (k*ds_bound_cylin(s)/(4i))*besselh(1,2,k*distance)*Jt(s)*(Xmid_bound_cylin(s)-x_far_field)/distance;
            else
                Hz_sc_far(I_far) = Hz_sc_far(I_far) - (k*ds_bound_cylin(s)/(4i))*besselh(1,2,k*distance)*Jt(s)*(Ymid_bound_cylin(s)-y_far_field)/distance;
            end
        end
        
    end
    
    for i = 1:Nx_domain % Calculating the Scattered Field and Total Field On the Entire Computational Doamin
        for j = 1:Ny_domain
            if ((X_domain(i)>Xmin_bound_cylin)&&(X_domain(i)<Xmax_bound_cylin)&&(Y_domain(j)>Ymin_bound_cylin)&&(Y_domain(j)<Ymax_bound_cylin))
                Hz_sc(i,j) = 0;
                Hz_inc(i,j) = 0;   
            else
                for s = 1:N_bound_cylin
                    distance = sqrt((X_domain(i)-Xmid_bound_cylin(s))^2 + (Y_domain(j)-Ymid_bound_cylin(s))^2);
                    if ((s>=1) && (s<=(Ny_bound_cylin-1))) % First segment, normal vector is x-hat
                        Hz_sc(i,j) = Hz_sc(i,j) - (k*ds_bound_cylin(s)/(4i))*besselh(1,2,k*distance)*Jt(s)*(X_domain(i)-Xmid_bound_cylin(s))/distance;
                    elseif ((s>=Ny_bound_cylin) && (s<=(Nx_bound_cylin+Ny_bound_cylin-2))) % Second segment, normal vector is y-hat
                        Hz_sc(i,j) = Hz_sc(i,j) - (k*ds_bound_cylin(s)/(4i))*besselh(1,2,k*distance)*Jt(s)*(Y_domain(j)-Ymid_bound_cylin(s))/distance;
                    elseif ((s>=(Nx_bound_cylin+Ny_bound_cylin-2)) && (s<=(Nx_bound_cylin+2*Ny_bound_cylin-3))) % Third segment, normal vector is -x-hat
                        Hz_sc(i,j) = Hz_sc(i,j) - (k*ds_bound_cylin(s)/(4i))*besselh(1,2,k*distance)*Jt(s)*(Xmid_bound_cylin(s)-X_domain(i))/distance;
                    else % Fourth segment, normal vector is -y-hat
                        Hz_sc(i,j) = Hz_sc(i,j) - (k*ds_bound_cylin(s)/(4i))*besselh(1,2,k*distance)*Jt(s)*(Ymid_bound_cylin(s)-Y_domain(j))/distance;
                    end
                end
                %Hz_inc(i,j) = exp(-Ima_Unit*k*X_domain(i)); % Incident Plane Wave Propagating in x-hat direction
                Hz_inc(i,j) = exp(-Ima_Unit*k*(X_domain(i)+Y_domain(j))/sqrt(2)); % Incident Plane Wave propagating in (x-hat + y-hat)/sqrt(2) direction
                Hz_total(i,j) = Hz_inc(i,j) + Hz_sc(i,j);
            end
        end
    end
    
    
    shift_count = 1; % Contour Angel Shifting Counter
    for n = 1:N_bound_cylin
        if (angle_bound_cylin(n)>min(angle_bound_cylin))
            shift_count = shift_count + 1; 
        else
            shift_count = shift_count - 1;
            break
        end
    end
    
    angle_bound_cylin = circshift(angle_bound_cylin,-1*shift_count);
    Xmid_bound_cylin = circshift(Xmid_bound_cylin,-1*shift_count);
    Ymid_bound_cylin = circshift(Ymid_bound_cylin,-1*shift_count);
    ds_bound_cylin = circshift(ds_bound_cylin,-1*shift_count);
    Jt = circshift(Jt,-1*shift_count);
    
    % Below Goes For Plotting Commands
    
    %plot(angle_bound_cylin,abs(Jt));
    %xlim([0 360]);
    %xlabel('Azimuthal Angle \phi   [\circ]','Fontsize',14);
    %ylabel('Induced Surface Current Density |J_t|    [A/m^2]','Fontsize',14);
    %title('Under Incident Plane Wave Propagating in $-\hat{x}$ (TE Polarization)  [ECE540 PJ3---Jianghuai Liu]','Interpreter','Latex','Fontsize',18);
    %title('Under Incident Plane Wave Propagating in $-\frac{\hat{x}+\hat{y}}{\sqrt{2}}$ (TE Polarization)  [ECE540 PJ3---Jianghuai Liu]','Interpreter','Latex','Fontsize',18);
    
    polarplot(angle_far_field,circshift(abs(Hz_sc_far),(N_angle_far_field-1)/2));
    %title('Far Field (at r=20$\lambda$) Scattered $|{H_z^{sc}}|$ Under Incident Plane Wave Propagating in $-\hat{x}$ (TE Polarization)  [ECE540 PJ3---Jianghuai Liu]','Interpreter','Latex','Fontsize',18);
    title('Far Field (at r=20$\lambda$) Scattered $|{H_z^{sc}}|$ Under Incident Plane Wave Propagating in $-\frac{\hat{x}+\hat{y}}{\sqrt{2}}$ (TE Polarization)  [ECE540 PJ3---Jianghuai Liu]','Interpreter','Latex','Fontsize',18);
    
    %[xPlot,yPlot] = meshgrid(X_domain./lambda,Y_domain./lambda);
    %contour(yPlot,xPlot,abs(Hz_total),6000);
    %xlabel('X/\lambda');
    %ylabel('Y/\lambda');
    %title('Total Magnetic Field $\left|H_z\right|$ in [A/m], with Incident Plane Wave Propagating in $-\hat{x}$ (TE Polarization)  [ECE540 PJ3---Jianghuai Liu]','Interpreter','Latex','Fontsize',14);
    %title('Total Magnetic Field $\left|H_z\right|$ in [A/m], with Incident Plane Wave Propagating in $-\frac{\hat{x}+\hat{y}}{\sqrt{2}}$ (TE Polarization)  [ECE540 PJ3---Jianghuai Liu]','Interpreter','Latex','Fontsize',14);
    %colorbar();
end

