function [A] = buildMatrix(panels)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
  N = length(panels);
  A = zeros(N);
  n = 0.5*ones(N,1);
  A(1:(N+1):end) = n;
   for i=1:max(N)
    for j=1:max(N)
      if (i ~= j)
        A(i,j) = (0.5/pi) * integratePanel(panels(i).x_center, ...
                                                 panels(i).y_center, ...
                                                 panels(j), ...
                                                 cos(panels(i).beta), ...
                                                 sin(panels(i).beta));
      end
    end
  end

end

