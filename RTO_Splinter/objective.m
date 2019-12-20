function f = eval_objective (x,J)
%OBJECTIVE Summary of this function goes here
  WG1 = x(1);
  WG2 = x(2);
  f = J.eval([WG1 WG2]);
end

