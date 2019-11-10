function [psi] = stream_function(strength, x_loc, y_loc, X, Y)
%stream_function(strength of stream, x location, y location, X,Y)

psi = strength ./ (2*pi) * arctan2((Y-y_loc),(X-x_loc));
end

