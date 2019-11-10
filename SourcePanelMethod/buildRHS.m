function [b] = buildRHS(panels, freestream)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
  
  N = length(panels);
  b = zeros(N);
  
  for i=1:(max(N))
    b(i) = -freestream.u_inf .* cos(freestream.alpha - panels(i).beta);
  end
  
end

