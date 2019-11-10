function [panel_object] = definePanel(x,y,N)
%definePanel - Discretize geometry using 'cosine' method


  r = (max(x) - min(x))/2; % radius of circle (chord of airfoil)
  x_center = (max(x) + min(x))/2; % x-coordinate center 
  
  % define x-coordinet of the circle points
  x_circle = x_center+r*cos(linspace(0.0,2*pi,N+1));
  
  x_ends = repmat(x_circle,1);
  y_ends = zeros(size(x_ends));
  
  I = 1;
  for i = 1:N
    while (I < (length(x)-1))
      if (((x(I)<=x_ends(i))&&(x_ends(i)<=x(I+1))) || ((x(I+1)<=x_ends(i))&&(x_ends(i)<= x(I))));
        break
      else
        I = I+1;
      end
    end
    a = (y(I+1) - y(I)) / (x(I+1) - x(I));
    b = y(I+1) - a * x(I+1);
    y_ends(i) = a * x_ends(i) + b;
    
  end
 
  for i = 1:max(N)
 
  panel_object(i) = createPanel(x_ends(i), y_ends(i), x_ends(i+1),y_ends(i+1));
  
  end
 
  
end

