% ME315 Numerical Problem
% Leslie T Rose
clear;clc;close all;

%% Setup Parameters

% Processor Properties

L = 75; %length in mm
W = 55; %Width in mm
H = 1;  %height in mm
N = 32; %number of cores
TDP = 250; %Total Dissipated Power
TdieMax = 68; %max allowed heat of die
f = 3; %processor frequency in GHz
fTurbo = 4.2; %processor turbo frequency in GHz
mem = 1;  %max Memory of system

k = 173; %W/mK thermal conductivity of CPU
kH2O = 0.61041; %W/mK Thermal Conductivity of Water
TDP = 250; %Total Dissipated Power

%Water inforation
Wmu = 8.5288e-4;
rho = 996.54;
R=W/L;

%Get nodal spacing, define dely, delV, and the arraySize
delx = input('Define nodal spacing (.5,1,2.5,5)mm: \n');
dely = delx*R;
delV = delx*dely;
meshSize = (round(L/delx))+1;


HeatGeneration; % run the Heat Generation Script
tic
Tamb=24;
Tdiemax = 68;

%options for Fsolve
options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-9,...
  'OptimalityTolerance',1e-9,'StepTolerance',1e-9,'MaxIter',1e5,...
  'MaxFunctionEvaluations',1e5);
%Setup anon function and then solve for Umin
Ueq = @(Uinf)TempSolver(Uinf,Tamb,delx,eGenAverage,meshSize);  %memoize function
memoizedUeq = memoize(Ueq);
Umin = fsolve(@(Uinf)memoizedUeq(Uinf)-Tdiemax,2.27,options); 


%create array of heat flux per square
q = dVarray.*(tempDist-Tamb).*repmat(h,meshSize,1);  
Q = sum(sum(q));
tempDistUmin = tempDist;

%Solve for 3*Umin
%Annon function to use with fsolve
%solve for Maximum ambient temperature 
Uinf=3*Umin; %Multiply Umin by 3 to model pump capacity
Tmaxeq = @(Tinf)TempSolver(Uinf,Tinf,delx,eGenAverage,meshSize); 
memTmaxeq = memoize(Tmaxeq);
TmaxAmb = fsolve(@(Tinf)memTmaxeq(Tinf)-.9*Tdiemax,35.57,options); 

%%

figure;
   hold on;
    Uminfig=contourf(xlimits,ylimits,tempDistUmin,100,'Linestyle', 'none');
    title('Temperature Distribution @ Minimum U')
    cBar = colorbar('southoutside'); cBar.Label.String = '°C';
    [y1,x1] = find(tempDistUmin==(max(max(tempDistUmin))))
    nx1=x1*75/meshSize;
    ny1=y1*(75/meshSize)*R;
    plot(nx1,ny1,'*');
  hold off;

figure;
  TempSolver(Uinf,TmaxAmb,delx,eGenAverage,meshSize);
  tempDistTambMax = tempDist;
  hold on;
    TambMaxCont=contourf(xlimits,ylimits,tempDist,100,'Linestyle','none');
    title('Temperature Distribution at Max Safe Ambient Temperature ');
    cBar = colorbar('southoutside'); cBar.Label.String = '°C';
    [y1,x1] = find(tempDistTambMax==(max(max(tempDistTambMax))))
    nx1=x1*75/meshSize;
    ny1=y1*(75/meshSize)*R;
    plot(nx1,ny1,'*');
  hold off;
  
%Solve for Reynolds number
ReyMin = (rho*Umin*L/1000)/Wmu
Rey3Min = (rho*Uinf*L/1000)/Wmu

%Sovle for Standard Deviation

contourDev = abs(tempDistTambMax-mean(mean(tempDistTambMax)));
figure('name','Standard Deviation Contour Map','Numbertitle','off');
 hold on;
  StdDev = contourf(xlimits,ylimits,contourDev,100,'Linestyle','none');
  title('Standard Deviation U = U_{min}');
  cBar = colorbar('southoutside'); cBar.Label.String = '°C';
 hold off;
 
MeanDev = mean(mean(abs(tempDistUmin-mean(mean(tempDistUmin)))));
MeanDev3U=mean(mean(abs(tempDistTambMax-mean(mean(tempDistTambMax)))));
%% Output Pertanent Information
%Clear Command Window and print out solutions 
clc;
convT=toc;

fprintf('Solutions have Converged in %6.2f seconds <-what do you think about them apples\n',convT);
fprintf('\nThe Optimal Flow velocity is %3.4f(m/s)',Umin);
fprintf('\nThe Max Safe Operational Ambient Temp is %2.3f°C\nTotal Q rejected is %3.2f(W)\n',TmaxAmb,Q);
fprintf('Standard Dev for 3U %3.2f °C\n',MeanDev3U);
fprintf('Reynolds @Umin = %3.3e \nReynolds @ 3*Umin = %3.3e\n',...
  ReyMin,Rey3Min);
fprintf('Standard Deviation is %1.3f °C\n',MeanDev);
