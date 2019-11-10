%% Code to visualize potential flows

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

%% Uniform flow 
U_inf = 1;
[u_free, v_free, psi_free] = freestream_functions(U_inf, N, Y);

% Create Source
source = source_flow(1,-1,0);
[u_source, v_source,psi_source] = source.functions(X,Y);

% Create sink
source2 = source_flow(-1,1,0);
[u_sink,v_sink,psi_sink] = source2.functions(X,Y);

% Create Doublet 
source1 = source_flow(10, -.25,0) ;
source2 = source_flow(-10, .25,0) ;
[u_source1,v_source1,psi_source1] = source1.functions(X,Y);
[u_sink1,v_sink1,psi_sink1] = source2.functions(X,Y);

%Create Vortex
vortex = vortexFlow(5,0,0);
[u_vortex, v_vortex,psi_vortex] = vortex.velocity(X, Y);


% Sink and source together
% Due to the equations being linear, superposition (addition), of the
% equations can be used to find the combined flow

fprintf(' 1 - Vortex, Source w/freestream \n');
fprintf(' 2 - Rankine Body \n');
fprintf(' 3 - Doublet \n');
fprintf(' 4 - Source w/freestream \n');
n = input('Enter number for type of potential flow: ');

switch n
  case 1 %vortex,source,freestream
      u =  u_source + u_vortex + u_free ;
      v =  v_source + v_vortex + v_free ;
      text = ('Circulation with Freestream');
  case 2 % Rankine body
      u = u_source + u_sink + u_free;
      v = v_source + v_sink + v_free;
      text = ('Rankine Body with freestream');
      
  case 3 % doublet
      u = u_source1 + u_sink1 + u_free;
      v = v_source1 + v_sink1 + v_free;
      text = (['Doublet with freestream']);
  
  case 4 % source only
      u = u_source + u_free;
      v = v_source + v_free;
      text = (['Source with freestream']);
    
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

