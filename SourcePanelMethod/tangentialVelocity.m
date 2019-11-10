function [panels] = tangentialVelocity(panels, freestream,sigma)

  N = length(panels);
  A2 = zeros(N);

  for i=1:max(N)
    for j=1:max(N)
      if (i ~= j)
        A2(i,j) = (0.5/pi) * integratePanel(panels(i).x_center, ...
                                            panels(i).y_center, ...
                                            panels(j), ...
                                           -sin(panels(i).beta), ...
                                            cos(panels(i).beta));
      end
    end
  end
  A2
  for i=1:(N)
    b2(i) = -freestream.u_inf * sin(panels(i).beta);
  end
  % A3 = A2*sigma;
  % tangential_velocity = A3+b2'
  tangential_velocity = A2*sigma+b2' %compute tangential velocity component
  
  for i=1:(max(N))
    panels(i).vt = tangential_velocity(i);  % add VT to objects
  end
end

