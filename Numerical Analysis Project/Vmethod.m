function [maxDie,A]= Vmethod(Uinf,Tamb,A,delx,eGenAverage)

% V method multi-mesh Gaussian 
% call the multimesh solver transpose and then repeat at the new meshsize

% V method  
n = [2 5 10 5 2]; % the number of iterations to be performed each coarseness

count = 0;
N = length(A);
A_Old = A;
while(1)
  
  A = multiMeshSolver(Uinf,Tamb,A,n(1),delx,eGenAverage);    % fine mesh
  
  A = resizeMesh(A);                    % resize down                     
  A = multiMeshSolver(Uinf,Tamb,A,n(2),delx,eGenAverage);    % medium mesh iteration
  
  A = resizeMesh(A);                    % resize down
  A = multiMeshSolver(Uinf,Tamb,A,n(3),delx,eGenAverage);    % coarse mesh iteration 
  
  A = repelem(A,2,2);                          % medium mesh transposition
  A = multiMeshSolver(Uinf,Tamb,A,n(4),delx,eGenAverage);    % medium mesh iteration
  
  A = repelem(A,2,2);                          % fine mesh transposition
  A = multiMeshSolver(Uinf,Tamb,A,n(5),delx,eGenAverage);    % fine mesh iteration

  error = abs((sum(sum(A)))-(sum(sum(A_Old))))/...
    (sum(sum(A_Old)));
  A_Old = A;
  if error < 1e-9  
  break
  end
  if mod(count,6.6e1)==0
    fprintf('.');
  end
  if mod(count, 2e2)==0
      clc;
      fprintf('Nodal Spacing of %1.1f mm, %5f Total Nodes\n Completed %6f iterations, Current error is %9e \n Processing',delx,N^2,count,error);
  end
    count = count+1;
end

maxDie = max(A);
end


