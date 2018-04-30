# This function takes an linear system of equations Ax=b with
# n being the dimension of the system, ep the maximum error of the expected 
# solutios, and itmax being the maximum number of iterations of the
# Gauss-Seidel algorithm.
#
# Warning: convergence is only guaranteed if the matrix A is diagonally dominant
# or if the matrix A is SPD

function x = gaussSeidel(A, b, n, ep, itmax)
  # Validar parametros de entrada
  if(n <= 0 || !ismatrix(A) || numel(A) != n*n || !iscolumn(b) || numel(b) != n)
    x = "ERR_INVALID_INPUT"; return; endif;
  # Tomamos x = [0, 0, 0, ..., 0] como o vetor inicial "x_0"
  x = xo = zeros(n, 1);
  for k = 1:itmax
    for i = 1:n
      x(1, i) += b(1, i) / A(i,i);
      for j = 1:(i-1)
        x(1, i) -= A(i, j) * x(1, j) / A(i, i);
      end for
      for j = (i+1):n
        x(1, i) -= A(i, j) * xo(1, j) / A(i, i);
      end for
    end for
    
    if(norm(x-xo, inf) < ep) return;
    xo = x;
  endfor
  x = "ERR_EXCEEDED_ITERATION";
endfunction

## Needs testing