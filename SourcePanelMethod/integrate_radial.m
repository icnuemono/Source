function [sum] = integrate_radial(p_i,p_j)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here

length = p_j.length;
int = @(s)(((p_i.x_center - (p_j.xa - sin(p_j.beta) .* s)) .* cos(p_i.beta) + ...
            (p_i.y_center - (p_j.ya + cos(p_j.beta) .* s)) .* sin(p_i.beta)) ./ ...
           ((p_i.x_center - (p_j.xa - sin(p_j.beta) .* s)).^2 + ...
            (p_i.y_center - (p_j.ya + cos(p_j.beta) .* s)).^2));
               
sum = integral(int,0.0,length);
               
end

