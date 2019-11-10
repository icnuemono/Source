function [u, v] = velocityField(panels, freestream, X, Y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

     u_free = freestream.u_inf * cos(freestream.alpha).* ones(size(X));
     v_free = freestream.u_inf * sin(freestream.alpha).* ones(size(X));
     
     N = length(panels);
     for i = 1:N
       for j = 1:N
       u_vel(i,j) = (panels(i).sigma / (2*pi)) * integratePanel(X(i),Y(j),panels(j),1,0);
       v_vel(i,j) = (panels(i).sigma / (2*pi)) * integratePanel(X(i),Y(j),panels(j),0,1);
       end
     end
     u_free
     v_free
     u_vel
     v_vel
     u = u_free + u_vel
     v = v_free + v_vel
end

