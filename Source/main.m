# This is the main script file of this assignment.
#
# It executes tests based on parameters received and computes the solution of
# the linear system using the Gauss-Seidel iterative method. 
format long

# First test (b))
disp("Primeiro Teste (b): ")
n = int32(50);
A = pentadiagMatrix(n);
b = sum(A)';
x = gaussSeidel(A, b, n, :, :)

# Second test (b))
disp("Segundo Teste (b): ")
n = int32(100);
A = pentadiagMatrix(n);
b = sum(A)';
x = gaussSeidel(A, b, n, :, :)

# Third test (c))
disp("Terceiro Teste (c): ")
n = int32(100);
A = pentadiagMatrix(n);
for(i=1:n)
  b(i) = 1/i;
endfor
b;
ep = 1e-10;
x = gaussSeidel(A, b, n, ep, :)