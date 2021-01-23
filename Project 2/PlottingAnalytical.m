
if (TE_Mode)
    Hz = zeros(Nnodes,1);
    H = Vmax(8);
    m = 3;
    n = 1;

    for i = 1:Nnodes
        Hz(i) = -H*cos(m*pi*x(i)/Length)*cos(n*pi*y(i)/Height);
    end

    scatter(x,y,10,Hz(:,1)-V(:,8));
    xlim([0 Length]);
    ylim([0 Height]);
    title('Difference Between Analytical and Numeric for H_{z}[A/m] in TE_{31}  ECE540,JH.L');
    colorbar();
    xlabel('X [m]');
    ylabel('X [m]');
end

if (TM_Mode)
    Ez = zeros(Nnodes,1);
    
    E = Vmax(7);
    m = 3;
    n = 2;

    for i = 1:Nnodes
        Ez(i) = E*sin(m*pi*x(i)/Length)*sin(n*pi*y(i)/Height);
    end

    scatter(x,y,10,Ez(:,1)-Vcomplete(:,7));
    xlim([0 Length]);
    ylim([0 Height]);
    title('Difference Between Analytical and Numeric for E_{z}[A/m] in TM_{32}  ECE540,JH.L');
    colorbar();
    xlabel('X [m]');
    ylabel('X [m]');
end

