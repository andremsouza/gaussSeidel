# This is the main script file of this assignment.
#
# It executes tests based on parameters received and computes the solution of
# the linear system using the Gauss-Seidel iterative method. 
format long

# First test (b))
disp("Primeiro Teste (b): \n")
n = uint32(50);
A = pentadiagMatrix(n);
b = sum(A)';
x = gaussSeidel(A, b, n, :, :, :)
disp("Primeiro Teste (b) (Erro): \n")
err = A*x - b

# Second test (b))
disp("Segundo Teste (b): \n")
n = uint32(100);
A = pentadiagMatrix(n);
b = sum(A)';
x = gaussSeidel(A, b, n, :, :, :)
disp("Segundo Teste (b) (Erro): \n")
err = A*x - b

# Third test (c))
disp("Terceiro Teste (c): \n")
n = uint32(100);
A = pentadiagMatrix(n);
for(i=1:n)
  b(i) = 1/i;
endfor
b;
ep = 1e-10;
x = gaussSeidel(A, b, n, :, ep, :)
disp("Terceiro Teste (c) (Erro): \n")
err = A*x - b