function [FCP_G,FCP_L] = Build_PVT_functions()
%BUILD_PVT_FUNCTIONS Summary of this function goes here
%   Detailed explanation goes here
load PVT.mat
load_libSplinter()

ZCG = reshape(CP_G',10,10)';
ZCL = reshape(CP_L',10,10)';


T = linspace(0,200,10)';
P = linspace(1,200,10)';

NP = length(P);
NT = length(T);

xs = zeros(NP*NT,2);
dCG = zeros(NP*NT,1);
dCL = zeros(NP*NT,1);
k = 1;
for i = 1:NP
   for j=1:NT
       xs(k,:) = [P(i) T(j)];
      dCG(k) = ZCG(i,j);
      dCL(k) = ZCL(i,j);
      k = k + 1;
   end
end

FCP_G = BSplineBuilder(xs, dCG, 4).build();
FCP_G.save('FCP_G.bspline')
FCP_L = BSplineBuilder(xs, dCL, 4).build();
FCP_L.save('FCP_G.bspline')
end

