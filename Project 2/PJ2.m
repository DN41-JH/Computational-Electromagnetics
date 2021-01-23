Length = 1; % Waveguide extentsion along x axis
Height = 0.5; % Waveguide extension along y axis
miu_r = 1; % Homogeneous magnetic permittivity of waveguide 
epsilon_r = 1; % Homogeneous electric permeability of waveguide

% Switches for choosing TE mode simulation or TM mode simulation
TE_Mode = 1; 
TM_Mode = 0;

% Call the mesh generation function to generate the triangular mesh, and
% obtain the connectivity array + nodes coordinate information
fd=@(p) drectangle(p,0,Length,0,Height);
[Coords,Connectivity]=distmesh2d(fd,@huniform,0.009,[0,0;Length,Height],[0,0;Length,0;0,Height;Length,Height]);

%fd=@(p) sqrt(sum(p.^2,2))-1;
%[Coords,Connectivity]=distmesh2d(fd,@huniform,0.03,[-1,-1;1,1],[]);

M = size(Connectivity,1); % Total Number of Elements
Nnodes = size(Coords,1); % Total Number of nodes

cTiny = 10^(-8);
NBoundaryNodes = 0;

% Calculate the number of boundary nodes
for n = 1:Nnodes
    if ((Coords(n,1)>(0-cTiny))&(Coords(n,1)<(0+cTiny)))
        NBoundaryNodes = NBoundaryNodes+1;
        continue
    elseif ((Coords(n,1)>(Length-cTiny))&(Coords(n,1)<(Length+cTiny)))
        NBoundaryNodes = NBoundaryNodes+1;
        continue
    elseif ((Coords(n,2)>(0-cTiny))&(Coords(n,2)<(0+cTiny)))
        NBoundaryNodes = NBoundaryNodes+1;
        continue
    elseif ((Coords(n,2)>(Height-cTiny))&(Coords(n,2)<(Height+cTiny)))
        NBoundaryNodes = NBoundaryNodes+1;
        continue
    end
end

iBoundaryNodes = zeros(NBoundaryNodes,1);
NInternalNodes = Nnodes - NBoundaryNodes; % Number of internal nodes
iInternalNodes = zeros(NInternalNodes,1);
i = 1;
j = 1;

% Store the Boundary node and Internal node index arrays
for n = 1:Nnodes
    if ((Coords(n,1)>(0-cTiny))&(Coords(n,1)<(0+cTiny)))
        iBoundaryNodes(i) = n;
        i = i+1;
        continue
    elseif ((Coords(n,1)>(Length-cTiny))&(Coords(n,1)<(Length+cTiny)))
        iBoundaryNodes(i) = n;
        i = i+1;        
        continue
    elseif ((Coords(n,2)>(0-cTiny))&(Coords(n,2)<(0+cTiny)))
        iBoundaryNodes(i) = n;
        i = i+1;        
        continue
    elseif ((Coords(n,2)>(Height-cTiny))&(Coords(n,2)<(Height+cTiny)))
        iBoundaryNodes(i) = n;
        i = i+1;       
        continue
    else
        iInternalNodes(j) = n;
        j = j+1;
    end
end

% Linear System: [K][Ez] = (k^2)[B][Ez] or [K][Hz] = (k^2)[B][Hz]
K_local = zeros(3,3,M); 
B_local = zeros(3,3,M);

K_global = zeros(Nnodes,Nnodes); 
B_global = zeros(Nnodes,Nnodes); 

ElementArea = zeros(1,M);
b = zeros(3,M);
c = zeros(3,M);

% Calculating the area of each triangular element
for e = 1:M
    b(1,e) = Coords(Connectivity(e,2),2) - Coords(Connectivity(e,3),2);
    b(2,e) = Coords(Connectivity(e,3),2) - Coords(Connectivity(e,1),2);
    b(3,e) = Coords(Connectivity(e,1),2) - Coords(Connectivity(e,2),2);
    c(1,e) = Coords(Connectivity(e,3),1) - Coords(Connectivity(e,2),1);
    c(2,e) = Coords(Connectivity(e,1),1) - Coords(Connectivity(e,3),1);
    c(3,e) = Coords(Connectivity(e,2),1) - Coords(Connectivity(e,1),1);
    ElementArea(e) = 0.5*(b(1,e)*c(2,e)-b(2,e)*c(1,e));
end

% Constructing each local [K] and [B] matrix
for e = 1:M
    B_local(:,:,e) = (epsilon_r*ElementArea(e)/12)*[2,1,1;1,2,1;1,1,2];
    for l = 1:3
        for k = 1:3
            K_local(l,k,e) = (1/(4*miu_r*ElementArea(e)))*(b(l,e)*b(k,e)+c(l,e)*c(k,e));
        end
    end
end

% Assembly process, adding the contribution of local [K] and [B] to the
% global [K] and [B]
for e = 1:M
    for l = 1:3
        for k = 1:3
            K_global(Connectivity(e,l),Connectivity(e,k)) = K_global(Connectivity(e,l),Connectivity(e,k))+K_local(l,k,e);
            B_global(Connectivity(e,l),Connectivity(e,k)) = B_global(Connectivity(e,l),Connectivity(e,k))+B_local(l,k,e);
        end
    end
end

% Plotting Preparation
x = zeros(Nnodes,1);
y = zeros(Nnodes,1);
for n = 1:Nnodes
    x(n) = Coords(n,1);
    y(n) = Coords(n,2);
end

if (TE_Mode)
    % Neumann Boundary Condition, No Treatment Needed
    % Diagnal of [D] is the eigenvalue
    % Each column of [V] is the field solution for Hz with each eigenvalue
    [V,D] = eig(K_global,B_global);
    D = sqrt(D);
    
    Vmax = zeros(Nnodes,1);
    for n = 1:Nnodes
        Vmax(n) = abs(max(V(:,n)));
    end
    
    
elseif (TM_Mode)
    AuxiliaryIBN = iBoundaryNodes;
    
    % Implementing Dirichlet Boundary Condition
    % Removing the rows and columns corresponding to the boundary nodes,
    % from the global [K] and [B] matrix
    for n = 1:NBoundaryNodes 
        K_global(AuxiliaryIBN(n),:) = [];
        K_global(:,AuxiliaryIBN(n)) = [];
        B_global(AuxiliaryIBN(n),:) = [];
        B_global(:,AuxiliaryIBN(n)) = [];
        AuxiliaryIBN = AuxiliaryIBN - ones(NBoundaryNodes,1);
    end
    
    % Diagnal of [D] is the eigenvalue
    % Each column of [V] is the field solution for Ez with each eigenvalue
    [V,D] = eig(K_global,B_global);
    D = sqrt(D);
    
    Vcomplete = zeros(Nnodes,Nnodes-NBoundaryNodes);
    TrackInternalNode = 1;
    TrackBoundaryNode = 1;
    
    % Complete the solved eigen-vector from (Nnodes-NBoundaryNodes) to
    % (Nnodes), by adding boundary nodes into it and assign field 0 on the
    % boundary nodes
    for i = 1:Nnodes
        if (iBoundaryNodes(TrackBoundaryNode)==i)
            Vcomplete(i,:) = 0;
            TrackBoundaryNode = TrackBoundaryNode + 1;
        else
            Vcomplete(i,:) = V(TrackInternalNode,:);
            TrackInternalNode = TrackInternalNode + 1;
        end
    end
       
    Vmax = zeros(Nnodes-NBoundaryNodes,1);
    for n = 1:(Nnodes-NBoundaryNodes)
        Vmax(n) = abs(max(Vcomplete(:,n)));
    end
end


