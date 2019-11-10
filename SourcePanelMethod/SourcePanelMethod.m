% Source-Panel Method on Arbitrary Airfoil
clear; clc;

% Import the airfoil from the .dat information *airfoil tools

airfoil_data = importdata('NACA0012.txt');  

for i = 1:size(airfoil_data)
    x_airfoil(i) = airfoil_data(i,1);
    y_airfoil(i) = airfoil_data(i,2);
end

%figure;
%plot(x,y)
%xlim([-0.1,1.1])
%ylim([-0.5,0.5])
%grid on
%pbaspect([1 .75 1])

% Discretize airfoil into panels

N = 50; 
panels = definePanel(x_airfoil,y_airfoil,N);

for i = 1:N
  x_ends(i) = panels(i).xa;
  y_ends(i) = panels(i).ya;
end


figure;
hold on
plot(x_airfoil,y_airfoil);
plot(x_ends,y_ends,'marker','o');
xlim([-0.1,1.1])
ylim([-0.5,0.5])
pbaspect([1 .75 1])
title('Airfoil with Discretized panels');
hold off

% Create freestream object
u_inf = 1.0;
alpha = 0.0;

freestream = Freestream(u_inf, alpha);

% Create the linear system
  A1 = buildMatrix(panels);
 b = buildRHS(panels, freestream);
  b1 = b(:,1);

% Solve system of linear equations 
sigma = linsolve(A1,b1);


for i=1:(N)
 % panels(i).vt = tangential_velocity(i);
  %panels(i).cp = 1.0 - (panels(i).vt / freestream.u_inf).^2;
  panels(i).sigma = sigma(i);  % add sigma value to objects
end

% Calculate Tangential Velocity
panels = tangentialVelocity(panels,freestream,sigma);

% Calculate Pressure Coefficient
panels = pressureCoefficient(panels,freestream);

figure;
hold on;
xc_u = [];
cp_u = [];
xc_l = [];
cp_l = [];

for i=1:N
  if (panels(i).loc == 'upper')
    xc_u(end+1) = panels(i).x_center;
    cp_u(end+1) = panels(i).cp;
  elseif (panels(i).loc == 'lower')
    xc_l(end+1) = panels(i).x_center;
    cp_l(end+1) = panels(i).cp;
  end
end
plot(xc_u,cp_u);
plot(xc_l,cp_l);
xlabel('x')
ylabel('Cp')
set(gca, 'YDir','reverse')
grid on
hold off;

figure;
hold on
plot(x_airfoil,y_airfoil);
plot(x_ends,y_ends,'marker','o');
xlim([-0.1,1.1])
ylim([-0.5,0.5])
pbaspect([1 .75 1])
title('Airfoil overlayed with Cp');
xc_u = [];
cp_u = [];
xc_l = [];
cp_l = [];

for i=1:N
  if (panels(i).loc == 'upper')
    xc_u(end+1) = panels(i).x_center;
    cp_u(end+1) = panels(i).cp;
  elseif (panels(i).loc == 'lower')
    xc_l(end+1) = panels(i).x_center;
    cp_l(end+1) = panels(i).cp;
  end
end
plot(xc_u,cp_u);
plot(xc_l,cp_l);
hold off

% Plot streamlines

% make meshgrid 


x_start= -1;
x_end = 2;
y_start = -0.3;
y_end = .3;

x = linspace(x_start, x_end,N);
y = linspace(y_start, y_end,N);
[X,Y] = meshgrid(x,y); 
Z = zeros(size(X));
w = ones(size(Y));

% calculate velocity field on the meshgrid

[u,v] = velocityField(panels, freestream, X, Y);


figure;
xstart = min(X)+(max(X)-min(X)).*rand(round(N/2),1); 
ystart = min(Y)+(max(Y)-min(Y)).*rand(round(N/2),1); 
zstart = min(Z)+(max(Z)-min(Z)).*rand(round(N/2),1);
hold on

xlim([x_start x_end]);
ylim([y_start y_end]);
h=streamslice(X,Y,u,v,3);
verts = stream2(X,Y,u,v,xstart,ystart);
set(h,'color','blue')
iverts = interpstreamspeed(X,Y,Z,u,v,w,verts,.025);
l = streamparticles(iverts,500,'animate',10,'ParticleAlignment','off','MarkerSize',3);
hold off

% figure;
% streamline(X,Y,u,v,xstart,ystart);

