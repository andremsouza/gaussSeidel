function A = pentadiagMatrix(n)
  A = diag(4*ones(1, n));
  for i = [-3,-1,1,3];
    A = A + diag(-ones(1, n - abs(i)), i);
  endfor
endfunction