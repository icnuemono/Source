function [u,v] = velocity_function(strength,x_loc, y_loc, X, Y)
%stream_function Creates potential flow (source, sink)
%   Provide x and y location for a sink use a negative strength

u = (strength / (2*pi)).*((X-x_loc)./((X-x_loc).^2+(Y-y_loc).^2));
v = (strength / (2*pi)).*((Y-y_loc)./((X-x_loc).^2+(Y-y_loc).^2));
end


