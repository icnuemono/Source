classdef Freestream
  %UNTITLED4 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    u_inf{mustBeNumeric}
    alpha{mustBeNumeric}
  end
  
  methods
    function obj = Freestream(u_inf, alpha)
      %UNTITLED4 Construct an instance of this class
      %   Detailed explanation goes here
      obj.u_inf = u_inf;
      obj.alpha = deg2rad(alpha);
    end

  end
end

