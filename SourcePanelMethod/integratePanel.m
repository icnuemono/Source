function [sum] = integratePanel(x,y,panels,dxdz,dydz)
%integratPanel determines strength at given panel center


length = panels.length;
int = @(s)(((x - (panels.xa - sin(panels.beta) .* s)) .* dxdz + ...
            (y - (panels.ya + cos(panels.beta) .* s)) .* dydz) ./ ...
           ((x - (panels.xa - sin(panels.beta) .* s)).^2 + ...
            (y - (panels.ya + cos(panels.beta) .* s)).^2));
           
sum = integral(int,0.0,length);
               
end

