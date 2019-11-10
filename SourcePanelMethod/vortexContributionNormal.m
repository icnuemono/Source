function [A] = vortexContributionNormal(panels,N)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


  A = zeros(N);
  %n = 0.5*ones(N_panels,1);
  %A(1:(N_panels+1):end) = n;
  for i=1:max(N)
    for j=1:max(N)
      if (i ~= j)
        A(i,j) = (0.5/pi) * integrate_tangential(panels(i).x_center, ...
                                              panels(i).y_center, ...
                                              panels(j), ...
                                              sin(panels(i).beta), ...
                                              -cos(panels(i).beta));
      end
    end
  end
end

