classdef vortexFlow
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    gamma
    x_pos
    y_pos
    X
    Y
    u
    v
    psi
    
  end
  
  methods
    function obj = vortexFlow(gamma,x_pos,y_pos)
      
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      obj.gamma = gamma;
      obj.x_pos = x_pos;
      obj.y_pos = y_pos;
      
    end
    
    function [u,v,psi] = velocity(obj,X,Y)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      u = obj.gamma / (2*pi) * (Y-obj.y_pos) ./ ((X-obj.x_pos).^2+(Y-obj.y_pos)^2);
      v = -obj.gamma / (2*pi) * (X-obj.x_pos)./ ((X-obj.x_pos).^2+(Y-obj.y_pos).^2);

      psi = obj.gamma ./ (4*pi) * log((X-obj.x_pos).^2+(Y-obj.y_pos)^2);
      
    end
  end
end

