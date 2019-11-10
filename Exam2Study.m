%% ME 482 Exam 2 Calculator

% NACA Code Solver

%NACA is XYZZ (page301)

%X = maximum camber in hundreths of chord
%Y = location of maximum camber in tenths
%ZZ = maximum thickness in hudreths of chord

c = 1.5;  % Chord length

X = 4; 
  camber = (X/100)*c
Y = 2;
  max_camber = (Y/10)*c
ZZ = 12;
  max_thickness = (ZZ/100)*c
  
 %%  Equations solve for based on variable 
 R = 287;
 cl = 1.2;
 T = 293;
 p = 101300;
 rho = p/(R*T)

syms Cl rho U_inf S

Lift = Cl*.5*rho*U_inf^2*S

L = solve(Lift,S)

%% Rankine Body

syms u v
