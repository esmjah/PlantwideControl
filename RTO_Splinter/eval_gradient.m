function g = eval_gradient (x,J)
  WG1 = x(1);
  WG2 = x(2);
  g = J.eval_jacobian([WG1 WG2]);
end
