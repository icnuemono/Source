% Leslie T Rose


fun = @rankineSolver; % create an anonymous function for 
x0 = [0.1,0.1];  % initial guesses for the solver 
x = fsolve(fun,x0)