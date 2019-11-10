classdef source_flow
  %UNTITLED13 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
  strength{mustBeNumeric}
  x_loc{mustBeNumeric}
  y_loc{mustBeNumeric}
  end
  
  methods
    function obj = source_flow(strength,x_loc, y_loc)
      if nargin == 3
        obj.strength = strength;
        obj.x_loc = x_loc;
        obj.y_loc = y_loc;
      end
    end
    
    function [u,v,phi] = functions(obj, X, Y)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      u = (obj.strength / (2*pi)).*((X-obj.x_loc)./((X-obj.x_loc).^2+(Y-obj.y_loc).^2));
      v = (obj.strength / (2*pi)).*((Y-obj.y_loc)./((X-obj.x_loc).^2+(Y-obj.y_loc).^2));
    
      phi = obj.strength ./ (2*pi) * atan2((Y-obj.y_loc),(X-obj.x_loc));
      
    end
    
    
  end
end

