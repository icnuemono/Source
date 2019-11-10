classdef createPanel
  %UNTITLED15 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    xa{mustBeNumeric}
    ya{mustBeNumeric} 
    xb{mustBeNumeric}
    yb{mustBeNumeric}
    sigma = 0;
    vt = 0;
    cp = 0;
    length;
    x_center;
    y_center;
    beta;
    loc;
    
  end
  
  methods
     function obj = createPanel(xa, ya, xb, yb)
      %UNTITLED15 Construct an instance of this class
      %   Detailed explanation goes here
       obj.xa = xa;
       obj.xb = xb;
       obj.ya = ya;
       obj.yb = yb;

        
       obj.x_center = (xa + xb)/2; % find center of x 
       obj.y_center = (ya + yb)/2;  %find center of y
       obj.length = sqrt((xb - xa)^2 + (yb - ya)^2); % find length of panel
       %Find orientation of  the panel
       
       if obj.xb - obj.xa <= 0
         obj.beta = acos((yb - ya)/obj.length);
       elseif  xb - xa > 0
         obj.beta = pi + acos(-(yb - ya) / obj.length);
       end
       
       % location of the panel with respect uppper or lower airfoil
       if (obj.beta <= pi)
         obj.loc = 'upper';
       else
         obj.loc = 'lower';
       end
       
       obj.sigma = 0.0;
       obj.vt = 0.0;
       obj.cp = 0.0;
       
    end


  end
end

