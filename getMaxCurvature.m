function [max_c] = getMaxCurvature (d, epsilon)
  max_c = 2*epsilon / (d*(1+epsilon))
endfunction
