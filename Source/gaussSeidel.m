# This function takes an linear system of equations Ax=b with
# n being the dimension of the system, ep the maximum error of the expected 
# solutios, and itmax being the maximum number of iterations of the
# Gauss-Seidel algorithm.
#
# Warning: convergence is only guaranteed if the matrix A is diagonally dominant
# or if the matrix A is SPD

function x = gaussSeidel(A=pentadiagMatrix(uint32(50)), b=(sum(A))', n=50, ep=1e-10, itmax=uint32(100))
  # Validar parametros de entrada
  if(n <= 0 || !ismatrix(A) || numel(A) != n*n || !iscolumn(b) || numel(b) != n || itmax <= 0 || ep <= 0)
    x = "ERR_INVALID_INPUT"; return; endif;
  # Tomamos x = [0, 0, 0, ..., 0] como o vetor inicial "x_0"
  x = xo = zeros(n, 1);
  for k = 1:itmax
    for i = 1:n
      x(i, 1) += b(i, 1) / A(i,i); # b_i term
      for j = 1:(i-1) # First summation
        x(i, 1) -= A(i, j) * x(j, 1) / A(i, i);
      endfor;
      for j = (i+1):n # Second summation
        x(i, 1) -= A(i, j) * xo(j, 1) / A(i, i);
      endfor;
    endfor;
    # If we have acceptable results, based on the ep variable
    if(norm(x-xo, inf) < ep) return; endif;
    xo = x;
  endfor;
  x = "ERR_EXCEEDED_ITERATION";
endfunction;

## Needs testing