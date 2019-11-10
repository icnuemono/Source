%% Freestream Solver for ME 482

clear; clc;

% Import the airfoil from the .dat information *airfoil tools

airfoil_data = importdata('NACA0012.txt');  

for i = 1:size(airfoil_data)
    x_airfoil(i) = airfoil_data(i,1);
    y_airfoil(i) = airfoil_data(i,2);
end

% %% Discretize airfoil into panels

N = 10; 
panels = definePanel(x_airfoil,y_airfoil,N);

for i = 1:N-1
  x_start(i) = panels(i).xa;
  y_start(i) = panels(i).ya;

end
 i = 0;
 x_end = [];
 y_end = [];
for i = N-2:N
  x_end(end+1) = panels(i).xb;
  y_end(end+1) = panels(i).yb;
end

figure;    % Plot original and discretized panels
hold on
plot(x_airfoil,y_airfoil)
plot(x_start,y_start,'marker','o')
plot(x_end,y_end,'marker','o')
xlim([-0.1,1.1])
ylim([-0.5,0.5])
pbaspect([1 .75 1])
hold off

% Create the meshgrid based on the x and y points on the discretized
% airfoil

x_start= -1;
x_end = 2;
y_start = -0.3;
y_end = .3;

% Grab the X and Y values of the points
for i = 1:N
  x1(i) = panels(i).x_center;
  y1(i) = panels(i).y_center;
end

x = linspace(x_start, x_end,N/2);
y = linspace(y_start, y_end,N/2);

x = [x x1];  %append x and y to include the points on the airfoil
y = [y y1];

x = sort(x);  %sort the values 
y = sort(y);

[X,Y] = meshgrid(x,y); %create a meshgrid using the points 

figure;    % Plot original and discretized panels
hold on
plot(x_airfoil,y_airfoil)
plot(X,Y,'b')
plot(x_start,y_start,'marker','o')
plot(x_end,y_end,'marker','o')
xlim([-0.1,1.1])
ylim([-0.5,0.5])
pbaspect([1 .75 1])
plot(X,Y)
hold off
% Create freestream object
u_inf = 30.0;
alpha = 0.0;

freestream = Freestream(u_inf, alpha);

% Calculate Flow 
  % Q = VA = freestream.alpha*A = psi
  
A = .4; %m
delx = .4/4;
dely = delx;
for i = 1:N
  if y(i)<=0
    d(i) = A+y(i);
    Q(i) = d(i)*freestream.u_inf;
  end
end
Q = flip(Q);
Q(end+1) = 0;
Q = flip(Q);
Q = [Q flip(Q)];
Q((round((length(y)+2)/2)))=[];
Q = Q';
Q = repmat(Q,1,(length(Q)));


solution = flowSolver(Q);

% I did not have time to go back and figure out the why my solutions were 
% not coming out as expected. I completed the problem an alternate method
% below

%% Array Method

%A matrix, created using a table
A = [-4,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;1,-4,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,1,-4,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,1,-4,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,1,-4,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,1,-4,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,1,-4,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;1,0,0,0,0,0,1,-4,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,1,-4,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,1,-4,1,0,0,0,0,0,1,0,0,0,0,0,0,0;0,0,0,1,0,0,0,0,0,1,-4,1,0,0,0,0,0,1,0,0,0,0,0,0;0,0,0,0,1,0,0,0,0,0,1,-4,1,0,0,0,0,0,1,0,0,0,0,0;0,0,0,0,0,1,0,0,0,0,0,1,-4,1,0,0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,1,0,0,0,0,0,1,-4,1,0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,0,1,0,0,0,0,0,1,-4,1,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,-4,1,0,0,0,0,0,1,0;0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,-4,1,0,0,0,0,0,1;0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,-4,1,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,-4,1,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,-4,1,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,-4,1,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,1,-4,1,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,1,-4,1;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,1,-4];
b = [-3;0;0;0;0;0;-3;6;0;0;0;0;0;-6;-9;-12;-24;-3;-6;-9;-9;12;12;12];

psi = inv(A)*b;  %Solve for psi

% Rearrange back into an array and then duplicate
psi = vec2mat(psi,7);
psi(end) = psi(4,3);
psi(4,3) = 12;
psi2 = flipud(psi);
psi2(1,:) = 12;
psi = vertcat(psi, psi2);
psi(4,:) = [];
psi(4,1:3) = 12;
psi(4,4:6) = 0;
psi(4,7) = 12  




% Still need to complete the u and v components 

x = [0:.10:.8];
y = [0:.10:.8];

% calculate inner matrix
for i = 2:6
  for j = 2:6
    u(i,j) = (psi(i,j+1) - psi(i,j-1))./2*dely;
    v(i,j) = -(psi(i+1,j) - psi(i-1,j))./2*delx;
  end
end

