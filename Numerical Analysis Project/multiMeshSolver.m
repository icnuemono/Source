function [A] = multiMeshSolver(Uinf,Tamb,A,i,delx,eGenAverage)

tempDist = A;  % Assign the tempDist as the passed in matrix. 
meshSize = length(A);
n = meshSize;
% CPU parameters
% These could be passed in, but this makes it easier. 
W = 55;               % width of the die
L = 75;               % length of the die
R = W/L;              % Ratio
k = 173;              % k w/m^2
H = 1;                % this should come from main
TdieMax = 68;         % Max desired temp of the Die in °C

%Properties of H20
mu = 8.5288e-4;       % Dynamic Viscocity Pa*s
rho = 996.54;         % density of water  kg/mm^3
kH2O = .61041;        % thermal conductivity of water w/mm*K
Cp = 4.1806*1000;     % thermal capacitance of water
x = 0:delx:L;         % length of die 75mm divided by delx

h = ones(1,meshSize);

for q = 1:meshSize
  Prandtl = mu*Cp/kH2O;
  Reynolds = rho.*Uinf.*(x(q)+125)/mu/1000;
  h(1,q) = 0.332.*kH2O.*power(Prandtl,1/3).*power(Reynolds,1/2)/((x(q)+125)/1000);
  
end
h = h./1e6;
assignin('base','h',h);
  %Solve for Bi

Bix = (h.*delx)/(k/1000);

%Setup Variables to be used in the equations.
eGenEq = (R^2*(delx^2)/(k/1000)).*eGenAverage;    % Setup eGen array

Tn = @(r,c,tempDist) tempDist((r-1),c);           % grab temp to the north
Ts = @(r,c,tempDist) tempDist((r+1),c);            % Grab temp to the south
Tw = @(r,c,tempDist) R^2*tempDist(r,(c-1));       % grab temp to the west
Te = @(r,c,tempDist) R^2*tempDist(r,(c+1));       % grab Temp to the east
Tinf = @(c) R^2*Bix(c)*delx*Tamb/H;               % Tinfinity Term
Denom = @(c) 2+2*R^2+(R^2*Bix(c)*delx/H);         % Denominator Value 


for i = 1:n              % So this will be run n number of times depending on nodal size 
  for c = 1:meshSize     
    for r = 1:meshSize
      
      % Solve for interior cells
      if r<meshSize && r>1 && c<meshSize && c>1
        tempDist(r,c) = (Tw(r,c,tempDist)+Te(r,c,tempDist)+Ts(r,c,tempDist)...
          +Tn(r,c,tempDist)+Tinf(c)+eGenEq(r,c))/Denom(c);
        %fprintf('Interior: %f\n', tempDist(r,c));
       
      % Solve left-edge
      elseif c==1 && r>1 && r<meshSize
        tempDist(r,c) = (Tinf(c)+Ts(r,c,tempDist)+Tn(r,c,tempDist)+...
          2*Te(r,c,tempDist)+eGenEq(r,c))/Denom(c); 
        %fprintf('Left Edge %f\n', tempDist(r,c));
        
      % TopEdge
      elseif c>1 && c<meshSize && r==1
        tempDist(r,c) = (Tw(r,c,tempDist)+Te(r,c,tempDist)+...
          2*Ts(r,c,tempDist)+Tinf(c)+eGenEq(r,c))/Denom(c);
        %fprintf('Top Edge %f\n', tempDist(r,c));
        
      % Bottom Edge
      elseif c>1 && c<meshSize && r==meshSize
        tempDist(r,c) = (Tw(r,c,tempDist)+Te(r,c,tempDist)+...
          2*Tn(r,c,tempDist)+Tinf(c)+eGenEq(r,c))/Denom(c);
        %fprintf('Bottom Edge %f\n', tempDist(r,c));  
        
      % Solve right-edge
      elseif c==meshSize && r>1 && r<meshSize
        tempDist(r,c) = (Tinf(c)+Ts(r,c,tempDist)+Tn(r,c,tempDist)+...
          2*Tw(r,c,tempDist)+eGenEq(r,c))/Denom(c);
        %fprintf('Right Edge %f\n', tempDist(r,c));
      
      % Top left Corner
      elseif r==1 && c==1
        tempDist(r,c) = (Tinf(c)+2*Ts(r,c,tempDist)+2*Te(r,c,tempDist)+...
          eGenEq(r,c))/Denom(c);
        %fprintf('Top Left %f\n', tempDist(r,c));
        
      % Bottom-Left Corner
      elseif c==1 && r==meshSize
        tempDist(r,c) = (Tinf(c)+2*Tn(r,c,tempDist)+2*Te(r,c,tempDist)+...
          eGenEq(r,c))/Denom(c);
        %fprintf('Bottom Left %f\n', tempDist(r,c));
      
      % Bottom-Right Corner
      elseif c==meshSize && r==meshSize
        tempDist(r,c) = (Tinf(c)+2*Tn(r,c,tempDist)+2*Tw(r,c,tempDist)+...
          eGenEq(r,c))/Denom(c); 
        %fprintf('Bottom Right %f\n', tempDist(r,c));
        
      % Top-Right Corner  
      elseif c==meshSize && r==1
        tempDist(r,c) = (Tinf(c)+2*Ts(r,c,tempDist)+2*Tw(r,c,tempDist)+...
          eGenEq(r,c))/Denom(c);
        %fprintf('Top Right %f\n', tempDist(r,c));
      end
    end
  end
  
A = tempDist;  
end