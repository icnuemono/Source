function [A] = sourceContributionNormal(panels,N)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here



  A = zeros(N);
  n = 0.5*ones(N,1);
  A(1:(N_panels+1):end) = n;
  for i=1:max(N_panels)
    for j=1:max(N_panels)
      if (i ~= j)
        A(i,j) = (0.5/pi) * integrate_radial(panels(i).x_center, ...
                                              panels(i).y_center, ...
                                              panels(j), ...
                                              cos(panels(i).beta), ...
                                              sin(panels(i).beta));
                                            
      end
    end
  end
end

