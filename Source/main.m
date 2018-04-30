# This is the main script file of this assignment.
#
# It receives user input from the terminal, and computes the solution of the
# linear system using the Gauss-Seidel iterative method. 
format long

while(1)
  try
    disp("Ax = b");
    # Receiving the dimension of the system and validating its value
    n = uint32(input("Digite a dimensao do sistema: "));
    ##if(n<=0) print("ERR_INVALID_INPUT\n\n"); continue; endif;
    
    # Receiving the matrix A, an custom user input, or the pentadiagMatrix(n).
    #
    # Warning: This script doesn't check if the matrix is diagonally dominant or
    # if the matrix is SPD.
    A = input("Digite a matriz A(ou ENTER, para usar a pentadiagonal): ");
    if(!numel(A)) A = pentadiagMatrix(n); endif;
    
    # Receiving vector b, and other parameters
    b = input("Digite o vetor coluna b do sistema: ")
    ep = input("Digite valor de epsilon (ε): ")
    itmax = uint32(input("Digite o numero maximo de iteraçoes do algoritmo: "))
    break;
  catch
    disp("ERR_INVALID_INPUT\n\n");
    continue;
  end_try_catch;
endwhile;

# Computing the solution, using the Gauss-Seidel iterative method
x = gaussSeidel(A, b, n, ep, itmax)