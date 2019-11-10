function [panels] = pressureCoefficient(panels, freestream)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
 N = length(panels);
 for i=1:(max(N))
  panels(i).cp = 1.0 - (panels(i).vt / freestream.u_inf).^2;
 end
end

