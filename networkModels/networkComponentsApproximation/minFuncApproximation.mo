within networkModels.networkComponentsApproximation;
function minFuncApproximation
input Real x1;
input Real x2;
output Real y;
algorithm
y := (-((x1-x2)^2+0.01^2)^0.5/2+(x1+x2)/2);
end minFuncApproximation;
