clear; clc; 

N = 50;

x_start= -5;
x_end = 5;
y_start = -2;
y_end = 2;

x = linspace(x_start, x_end,N);
y = linspace(y_start, y_end,N);
[X,Y] = meshgrid(x,y); 
Z = zeros(size(X));
w = ones(size(Y));

%  Create line of sources 

posx = linspace(x_start,x_end,5);
posy = [0 0 0 0 0];
u = zeros(N);
v = zeros(N);

for i = 1:5
  source(i) = source_flow(1, posx(i),posy(i));
  [u_source(i),v_source(i)] = source(i).functions(X,Y);
  u = u + source(i);
  v = v + source(i);
end

figure;
xstart = min(X)+(max(X)-min(X)).*rand(round(N/2),1); 
ystart = min(Y)+(max(Y)-min(Y)).*rand(round(N/2),1); 
zstart = min(Z)+(max(Z)-min(Z)).*rand(round(N/2),1);
hold on
xlim([x_start x_end]);
ylim([y_start y_end]);
title(text);
h=streamslice(X,Y,u,v,3);
verts = stream2(X,Y,u,v,xstart,ystart);
xstag = source.strength/(2*pi*U_inf);
x2stag = source2.strength/(2*pi*U_inf);
scatter(xstag,0,20,'red') 
scatter(x2stag,0,20,'red')
set(h,'color','blue');
iverts = interpstreamspeed(X,Y,Z,u,v,w,verts,.025);
streamparticles(iverts,500,'animate',10,'ParticleAlignment','off','MarkerSize',3);
hold off

  
  