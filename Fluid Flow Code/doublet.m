function [u,v, psi] = doublet(strength,x_source, y_source, X,Y)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
u = (-strength / (2*pi) * ((X-x_source).^2 - (Y-y_source).^2) / ((X-x_source).^2 + (Y-y_source).^2).^2);
v = (-strength / (2*pi) * 2 * (X-x_source)*(Y-y_source) / ((X-x_source).^2+(Y-y_source).^2)^2);
end

