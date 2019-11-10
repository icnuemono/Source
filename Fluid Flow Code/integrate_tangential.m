function [sum] = integrate_tangential(p_a,p_b)
%UNTITLED18 Summary of this function goes here
%   Detailed explanation goes here
length = p_b.length;
int = @(s)((-(p_a.x_center - (p_b.xa - sin(p_b.beta) .* s)) .* sin(p_a.beta) + ...
                  (p_a.y_center - (p_b.ya + cos(p_b.beta) .* s)) .* cos(p_a.beta)) ./ ...
                 ((p_a.x_center - (p_b.xa - sin(p_b.beta) .* s)).^2 + ...
                 (p_a.y_center - (p_b.ya + cos(p_b.beta) .* s)).^2));
               
sum = integral(int,0.0,length);
end

