function [S] = pointIterative(x)

S(1) = 2*x(1)+x(2)+x(3)-5;
S(2) = -x(1)+3*x(2)-x(3)-2;
S(3) = x(1) - x(2) + 2*x(3) -5;

end