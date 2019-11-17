within networkModels.networkComponents;
function heaviside
  input Real x;
  output Real y;
algorithm
  y := 1 / (1 + exp(-2 * 1 * x));
end heaviside;
