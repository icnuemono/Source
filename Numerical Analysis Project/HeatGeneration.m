%% Heat Generation Array 
%so the eGen function will stay the same muther fudger... What needs to
%change is the limits of integration in the eGenInt variable... 

eGenAv = zeros(meshSize);
xlimits = 0:delx:75;
ylimits = 0:dely:55;


%start computing 
fprintf('Beginning Heat Generation Calculation\n');

tic %stat timer
egenPerA=@(x,y)(integral2(@(x,y)tempCalc(x,y),0,L,0,W,'AbsTol',0,...
  'RelTol',1e-25));
A=TDP/egenPerA(L,W);

%%
eGenFun = @(x,y)A*tempCalc(x,y);
dumbx = 0:0.1:L;
dumby = 0:0.1:W;
[xx,yy] = meshgrid(dumbx,dumby);
figure;
%Output the Volumetric Heat Generation Contour MAp
contourf(xx,yy, eGenFun(xx,yy),40,'Linestyle','none')                       
title('Volumetric Heat Generation')
cBar = colorbar('southoutside'); 
cBar.Label.String = 'Color bar label ?units';

for c=1:(meshSize)  
  for r=1:(meshSize)     
     if c==1
       alpha=0;beta=1;
     elseif c==meshSize
       alpha = 1; beta = 0;
     else
       alpha = 1; beta = 1;
     end
     if r==1
       gamma = 0; delta = 1;
     elseif r==meshSize
       gamma = 1;delta = 0;
     else
       gamma = 1; delta = 1; 
     end     
     x0 = xlimits(c)-alpha*(delx/2);
     x1 = xlimits(c)+beta*(delx/2);
     y0 = ylimits(r)-gamma*(dely/2);
     y1 = ylimits(r)+delta*(dely/2);
   
     eGenAv(r,c)=integral2(@(x,y)A*tempCalc(x,y),x0,x1,y0,y1,'AbsTol',0,...
       'RelTol',1e-25);
  end   
end

%now do the division by volume 
dVarray = delx*dely*ones(meshSize);
dVarray(1,:)=dVarray(1,:)/2;
dVarray(meshSize,:)=dVarray(meshSize,:)/2;
dVarray(:,1)=dVarray(:,1)/2;
dVarray(:,meshSize)=dVarray(:,meshSize)/2;

eGenAverage=eGenAv./dVarray;
%eGenRend = contourf(xlimits,ylimits,eGenAv,40,'Linestyle','none');
fprintf('Heat Generation Array Complete\n');
toc %end timer
%


