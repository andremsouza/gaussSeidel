# This function builds a pentadiagonal matrix with dimensions n x n
function A = pentadiagMatrix(n) # n > 0; n is an int32 variable
  if(!isinteger(n) || n <= 0) # check invalid inputs and return error
    A = "ERR_INVALID_INPUT"; return; endif;
  A = diag(4*ones(n, 1)); # main diagonal matrix
  for i = [-3,-1,1,3] # computing the "shifted" diagonals
    if(n > abs(i))
      A = A + diag(-ones(1, n - abs(i)), i); endif;
  endfor;
endfunction;