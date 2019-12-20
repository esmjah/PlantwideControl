function [WL,WG,PWH,TWH] = Build_well_functions(Pres,fineGrid)
%BUILD_WELL160_FUNCTIONS Summary of this function goes here
%   Detailed explanation goes here

if(fineGrid==1)
    filename = ['Well_' num2str(Pres) '_Olga2014_fine_Data.mat'];
else
     filename = ['Well_' num2str(Pres) '_Olga2014_Data.mat'];
end

load(filename);

NP = length(Pwh);
NW = length(WGin);

xs = zeros(NP*NW, 2);
P_whds = zeros(NP*NW, 1);
T_whds = zeros(NP*NW, 1);
WG_whds = zeros(NP*NW, 1);
WL_whds = zeros(NP*NW, 1);
k = 1;
for i = 1:NP
   for j = 1:NW
      xs(k,:) = [Pwh(i) WGin(j)];
      P_whds(k) = P_whd(i,j);
      T_whds(k) = T_whd(i,j);
      WG_whds(k) = WG_whd(i,j);
      WL_whds(k) = WL_whd(i,j);
      k = k + 1;
   end
end

WL = BSplineBuilder(xs, WL_whds, 4).build();

WG = BSplineBuilder(xs, WG_whds, 4).build();

PWH = BSplineBuilder(xs, P_whds, 4).build();

TWH = BSplineBuilder(xs, T_whds, 4).build();

end

