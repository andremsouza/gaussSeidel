# This function takes an linear system of equations Ax=b with
# n being the dimension of the system, ep the maximum error of the expected 
# solutios, and itmax being the maximum number of iterations of the
# Gauss-Seidel algorithm.
#
# Warning: convergence is only guaranteed if the matrix A is diagonally dominant
# or if the matrix A is SPD

function x = gaussSeidel(A, b, n=3, xo = zeros(n, 1), ep=1e-10, ...
                          itmax=uint32(1048576))
  # Validate all input parameters
  if(n <= 0 || !ismatrix(A) || numel(A) != n*n || !iscolumn(b) || ...
      numel(b) != n || itmax <= 0 || ep <= 0 || !iscolumn(xo) || numel(xo) != n)
    x = "ERR_INVALID_INPUT"; return; endif;
  x = xo;
  for k = 1:itmax
    xo = x;
    for i = 1:n
      j = 1:n; # To use in the summation
      j(i) = []; # Removing j = i
      xAux = x; # To use in the summation
      xAux(i) = []; # Removing j = i
      x(i) = (b(i) - sum(A(i, j) * xAux)  )/A(i, i);
    endfor;
    # If we have acceptable results, based on the ep variable
    err = norm(x-xo, inf);
    if(err < ep) return; endif;
  endfor;
  # If it gets at this point, the function has exceeded the itmax value
  x = "ERR_EXCEEDED_ITERATION";
endfunction;