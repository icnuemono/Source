%Gauss-seidel iteration solver for flow problem 
function [Psi]= flowSolver(Q);

N = length(Q)


Psi_Old=Q;
Psi=Q;


count = 0;

for i = 1:1e6
  for c = 1:N     
    for r = 1:N
      
      % Solve for interior cells
      if r<N && r>1 && c<N && c>1
        Psi(r,c) = (Psi(r,c+1)+Psi(r,c-1)+Psi(r+1,c)...
          +Psi(r-1,c)-4*Psi(r,c));
        
      end
    end
   
  end
 
  error = abs((sum(sum(Psi)))-(sum(sum(Psi_Old))))/...
    (sum(sum(Psi_Old)));
  Psi_Old = Psi;
  if error < 1e-9  
    break
  end
  if mod(count,6.6e1)==0
    fprintf('.');
  end
  if mod(count, 2e2)==0
      clc;
      fprintf(' %5f Total Nodes\n Completed %6f interations, Current error is %9e \n Processing',N^2,count,error);
  end
    count = count+1;
end


end
