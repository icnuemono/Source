%% Discretizing of A circle
clear; clc;
tic;
u_inf = 1 ;  %freestream velocity
r = 1.0;  %radius of the circle
[x_center, y_center] = deal(0,0);
theta = linspace(0.0,2*pi,100);
[x_cylinder,y_cylinder] = deal(x_center+r*cos(theta),y_center+r*sin(theta));

N_panels = 10;


%find the end points of the panels

x_ends = r * cos(linspace(0.0,2*pi,N_panels +1 ));
y_ends  = r * sin(linspace(0.0,2*pi,N_panels + 1));

% define the panels



for i = 1:max(N_panels)
 
  panels(i) = Panels(x_ends(i), y_ends(i), x_ends(i+1),y_ends(i+1));
  
end

figure;
hold on;
grid on;
xlim([-r-0.2,r+0.2]);
ylim([-r-0.2,r+0.2]);
plot(x_cylinder,y_cylinder);
plot(x_ends, y_ends);
for i = 1:max(N_panels)
  scatter(panels(i).x_center,panels(i).y_center)
end

hold off;
pbaspect([1 1 1])


%% Flow-tangency bounday condition

A = zeros(N_panels);
n = 0.5*ones(N_panels,1);
A(1:(N_panels+1):end) = n;
for i=1:max(N_panels)
  for j=1:max(N_panels)
    if (i ~= j)
      A(i,j) = (0.5/pi) * integrate_radial(panels(i),panels(j));
    end
  end
end

for i=1:(max(N_panels))
  b(i) = -u_inf .* cos(panels(i).beta);
end
 b = b';


sigma = linsolve(A,b);

%% Solving for Pressure Coefficients 


A2 = zeros(N_panels);
%n = 0.5*ones(N_panels,1);
%A(1:(N_panels+1):end) = n;
for i=1:max(N_panels)
  for j=1:max(N_panels)
    if (i ~= j)
      A2(i,j) = (0.5/pi) * integrate_tangential(panels(i),panels(j));
    end
  end
end

for i=1:(max(N_panels))
  b2(i) = -u_inf * sin(panels(i).beta);
end

b2 = b2';
tangential_velocity = (A2*sigma)+b2;

for i=1:(max(N_panels))
  panels(i).vt = tangential_velocity(i);
  panels(i).cp = 1.0 - (panels(i).vt / u_inf).^2;
  panels(i).sigma = sigma(i);
end



cp_analytical = 1.0-4.*(y_cylinder./r).^2;
 
figure;
hold on
plot(x_cylinder, cp_analytical)
for i=1:(max(N_panels))
scatter(panels(i).x_center, panels(i).cp)
end
hold off
toc;
  
  
  