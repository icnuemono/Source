%% Nozzle Design for Liquid Fueled Rocket
clear;clc;


F = 10;     % desired thrust
Isp = 244;  % specific Impulse
Pc = 200;   % chamber pressure in PSI (desired)
r = 1;      % Fuel/oxygen mixture ratio
wt = F/Isp; % fuel flow 
wo = wt*r/(r+1); % flow rate of the oxydizer
wf = wt*r/(r+1); % flow rate of the fuel 


% Parameters

M = 24 ;         % Molecular weight of the gas
R_bar = 1545.32; % ft-lb/lb(deg)R
R = R_bar/M ;    % gas constant
gc = 32.2;       % ft/sec^2
gamma = 1.2;     %
Pe = 100;        % pressure in the combustion chamber 
Patm = 14.7;     % local atmospheric pressure (could end up being formula to find for altitude)
Tf = 5180;       % combustion flame temperature in F�
Te = 25;         % Temperature and nozzle exit 

% The goal is accelerating the gas to the local speed of sound M=1

Tc = Tf+460; % combustion flame temperature in degree Rankine
Tt = Tc*(1/(1+((gamma-1)/2))); % Temperature of gases at nozzle throat.
Pt = Pe*(1+((gamma-1)/2)); % gas pressure at the throat nozzle 

Me = (2/(gamma-1))*(((Pe/Patm)^((gamma-1)/gamma))-1) % mach at nozzle exit 
Me = sqrt(Me)
At = wt/(Pt*sqrt((R*Tt)/(gamma*gc))) % Throat cross-secitonal area

Ae = (At/Me)*((1+((gamma-1)/2)*Me^2)/((gamma+1)/2))^((gamma+1)/(2*(gamma-1))) % Nozzle exit area

A_ratio = At*(Ae/At)
T_ratio = Tc*(Te/Tc)