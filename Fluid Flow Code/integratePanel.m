function [sum] = integratePanel(x,y,panel,dxdz,dydz)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here

length = panel.length;
int = @(s)(((x - (panel.xa - sin(panel.beta) .* s)) .* dxdz + ...
            (y - (panel.ya + cos(panel.beta) .* s)) .* dydz) ./ ...
           ((x - (panel.xa - sin(panel.beta) .* s)).^2 + ...
            (y - (panel.ya + cos(panel.beta) .* s)).^2));
               
sum = integral(int,0.0,length);
               
end

