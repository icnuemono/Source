% Leslie T Rose
% Rankine Body Solver for Exam

fun = @rankineSolver; % create an anonymous function for 
x0 = [0.1,0.1];  % initial guesses for the solver 
x = fsolve(fun,x0)  % calll fsolve and solve for Gamma and a

function [F] = rankineSolver(x)
%rankineSolver(x) used to find gamma and a of rankine body
% x(1) = a  x(2) = gamma 
h = 0.03;
l = 0.05;
Vinf = 22;
F(1) = (x(1) * cot((Vinf*h*pi) / x(2)))-h;
F(2) = (x(1) * sqrt((x(2)/(x(1)*pi*Vinf))+1))-l;

end

