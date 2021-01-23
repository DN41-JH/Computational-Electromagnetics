epsilon = 8.8542*10^(-12);
miu = 4*pi*10^(-7);
v = sqrt(1/(miu*epsilon));
impedance = sqrt(miu/epsilon);

Wavelength = 0.01; % Wavelength of the single-frequency current source
Frequency = v/Wavelength; % Corresponding frequency of the single-frequency current source
k = 2*pi/Wavelength; % Wave number
tao_p = 3*(1/Frequency);
X = 5*Wavelength;
Y = 5*Wavelength;
X_c = 0.5*X;
Y_c = 0.5*Y;
dX = Wavelength/40;
dY = Wavelength/40;
nX = (X/dX)+1;
nY = (Y/dY)+1;
dT = 1/(v*sqrt((1/(dX)^2)+(1/(dY)^2)));
L_X = 40;
L_Y = 40;
TurnOnPML = 1;
nStep = 800;

Bx = zeros(nStep,nX,nY+1);
By = zeros(nStep,nX+1,nY);
Hx = zeros(nStep,nX,nY+1);
Hy = zeros(nStep,nX+1,nY);
Dz = zeros(nStep,nX,nY);
Ez = zeros(nStep,nX,nY);
Ez_analytical = zeros(nX,nY);
Jz = zeros(nStep,nX,nY);

m = 2;
SigmaX = zeros(nX,nY);
SigmaY = zeros(nX,nY);
SigmaZ = 0;
SigmaX_max = -((m+1)/(2*impedance*L_X*dX))*log(0.01);
SigmaY_max = -((m+1)/(2*impedance*L_Y*dY))*log(0.01);

if (TurnOnPML) 
    for i = 1:L_X+1
        for j = 1:nY
            SigmaX(i,j) = SigmaX_max*(abs(i-(L_X+1))/L_X)^m;
        end
    end
    for i = (nX-L_X):nX
        for j = 1:nY
            SigmaX(i,j) = SigmaX_max*(abs(i-(nX-L_X))/L_X)^m;
        end
    end
    for i = 1:nX
        for j = 1:L_Y+1
            SigmaY(i,j) = SigmaY_max*(abs(j-(L_Y+1))/L_Y)^m;
        end
    end
    for i = 1:nX
        for j = (nY-L_Y):nY
            SigmaY(i,j) = SigmaY_max*(abs(j-(nY-L_Y))/L_Y)^m;
        end
    end
end


AlphaX = zeros(nX,nY);
AlphaY = zeros(nX,nY);
BetaX = zeros(nX,nY);
BetaY = zeros(nX,nY);
AlphaZ = (epsilon/dT)-(SigmaZ/2);
BetaZ = (epsilon/dT)+(SigmaZ/2);

for i = 1:nX
    for j = 1:nY
        AlphaX(i,j) = (epsilon/dT)-(SigmaX(i,j)/2);
        BetaX(i,j) = (epsilon/dT)+(SigmaX(i,j)/2);
        AlphaY(i,j) = (epsilon/dT)-(SigmaY(i,j)/2);
        BetaY(i,j) = (epsilon/dT)+(SigmaY(i,j)/2);
    end
end

for n = 1:nStep
   Jz(n,(nX-1)/2,(nY-1)/2) = (1-exp(-((n-0.5)*dT/tao_p)))*sin(2*pi*Frequency*(n-0.5)*dT);
   %Jz(n,(nX-1)/2,(nY-1)/2) = sin(2*pi*Frequency*(n-0.5)*dT);
   %Jz(n,(nX-1)/2,60) = (1-exp(-((n-0.5)*dT/tao_p)))*sin(2*pi*Frequency*(n-0.5)*dT);
end

%plot(Jz(:,(nX-1)/2,(nY-1)/2));

for n = 2:nStep
    for i = 1:nX
        for j = 2:nY
            Bx(n,i,j) = (1/((BetaY(i,j-1)+BetaY(i,j))/2))*((((AlphaY(i,j-1)+AlphaY(i,j))/2))*Bx(n-1,i,j)-(epsilon/dY)*(Ez(n-1,i,j)-Ez(n-1,i,j-1)));
        end
    end
    
    for i = 2:nX
        for j = 1:nY
            By(n,i,j) = (1/BetaZ)*(AlphaZ*By(n-1,i,j)+(epsilon/dX)*(Ez(n-1,i,j)-Ez(n-1,i-1,j)));
        end
    end
    
    for i = 1:nX
        for j = 2:nY
            Hx(n,i,j) = (1/BetaZ)*(AlphaZ*Hx(n-1,i,j)+(((BetaX(i,j-1)+BetaX(i,j))/2)/miu)*Bx(n,i,j)-(((AlphaX(i,j-1)+AlphaX(i,j))/2)/miu)*Bx(n-1,i,j));
        end
    end
    
    for i = 2:nX
        for j = 1:nY
            Hy(n,i,j) = (1/((BetaX(i-1,j)+BetaX(i,j))/2))*(((AlphaX(i-1,j)+AlphaX(i,j))/2)*Hy(n-1,i,j)+(((BetaY(i-1,j)+BetaY(i,j))/2)/miu)*By(n,i,j)-(((AlphaY(i-1,j)+AlphaY(i,j))/2)/miu)*By(n-1,i,j));
        end 
    end
    
    for i = 2:nX-1
        for j = 2:nY-1
            Dz(n,i,j) = (1/BetaX(i,j))*(AlphaX(i,j)*Dz(n-1,i,j)+(epsilon/dX)*(Hy(n,i+1,j)-Hy(n,i,j))-(epsilon/dY)*(Hx(n,i,j+1)-Hx(n,i,j))-epsilon*Jz(n,i,j));
        end
    end
    
    for i = 2:nX-1
        for j = 2:nY-1
            Ez(n,i,j) = (1/BetaY(i,j))*(AlphaY(i,j)*Ez(n-1,i,j)+(1/epsilon)*(BetaZ*Dz(n,i,j)-AlphaZ*Dz(n-1,i,j)));
        end
    end
    
    %Ez(n,1:98,140) = 0;  %For Single Slit
    %Ez(n,102:nX,140) = 0; % For Single Slit
    
    Ez(n,1:60,120) = 0; % For Double Slit
    Ez(n,64:130,120) = 0; % For Double Slit
    Ez(n,134:nX,120) = 0; % For Double Slit
end

%for i = 1:nX
%    for j = 1:nY
%        Ez_analytical(i,j) = 10^(-7.3)*real(-((k^2)/(4*2*pi*Frequency*epsilon))*besselh(0,2,k*abs(sqrt((((i-1)*dX)^2)+(((j-1)*dY)^2))-sqrt(((0.5*X)^2)+((0.5*Y)^2))))*exp(-1i*2*pi*Frequency*dT*50));
%    end
%end

[xPlot,yPlot] = meshgrid(0:dX:X,0:dY:Y);
[xDomain,yDomain] = meshgrid(dX*(L_X+1):dX:dX*(nX-L_X),dY*(L_Y+1):dY:dY*(nY-L_Y));

%contour(xPlot,yPlot,reshape(Ez(150,:,:),[nX,nY]),2000);
%contour(xPlot,yPlot,reshape(Hx(1000,1:nX,1:nY),[nX,nY]),2000);
%contour(xPlot,yPlot,reshape(Hy(1000,1:nX,1:nY),[nX,nY]),2000);
%contour(xDomain,yDomain,reshape(Ez(1000,L_X+1:nX-L_X,L_Y+1:nY-L_Y),[nX-2*L_X,nY-2*L_Y]),2000);
