function [B] = resizeMesh(A)

B = ones(length(A)/2);
N = length(B);

for i = 1:N
  for j = 1:N
    
      B(i,j) = (A(i,j)+A(i+1,j)+A(i,j+1)+A(i+1,j+1))/4;
    
  end
end
