function [u, v, psi] = vortex_function(strength,x_pos,y_pos,X,Y)
%vortex_functioncalculates the u,v, and psi components of a vortex
%   Detailed explanation goes here

u = strength / (2*pi) * (Y-y_pos) ./ ((X-x_pos).^2+(Y-y_pos)^2);
v = -strength / (2*pi) * (X-x_pos)./ ((X-x_pos).^2+(Y-y_pos).^2);

psi = strength ./ (4*pi) * log((X-x_pos).^2+(Y-y_pos)^2);
end

