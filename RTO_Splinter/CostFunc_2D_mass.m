function J = CostFunc_2D_mass(x,Lambda,WL_w1,WL_w2)
%COSTFUNC Summary of this function goes here
%   Detailed explanation goes here

Winj1 = x(1);
Winj2 = x(2);
Delta_Pwh1 = x(3);
Delta_Pwh2 = x(4);
Pm = x(6);

Pwh1 = Pm + Delta_Pwh1;
Pwh2 = Pm + Delta_Pwh2;

PriceL = Lambda(1);
PriceG = Lambda(2);

WL1 = fnval(WL_w1,{Pwh1,Winj1});


WL2 = fnval(WL_w2,{Pwh2,Winj2});

J = (-PriceL*(WL1+WL2) + PriceG*(Winj1+Winj2));




end

