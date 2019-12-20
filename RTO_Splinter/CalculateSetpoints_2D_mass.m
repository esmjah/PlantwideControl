function [WLtot,alpha1,alpha2,alpha,P1,P2,Pin,T1,T2,T_INL,T_top] = CalculateSetpoints_2D_mass( x,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L,FH_G,FH_L,FTH_TOP)
%NonlinCon Summary of this function goes here
%   Detailed explanation goes here




WG_w1 = w1.FWG;
WG_w2 = w2.FWG;

WL_w1 = w1.FWL;
WL_w2 = w2.FWL;

TWH_w1 = w1.FTWH;
TWH_w2 = w2.FTWH;

PWH_w1 = w1.FPWH;
PWH_w2 = w2.FPWH;

Winj1 = x(1);
Winj2 = x(2);
Delta_Pwh1 = x(3);
Delta_Pwh2 = x(4);
Delta_Ptop = x(5);

Pwh1 = x(6) + Delta_Pwh1;
Pwh2 = x(6) + Delta_Pwh2;
Ptop = 5e5 + Delta_Ptop;


T1 = fnval(TWH_w1,{Pwh1,Winj1});
T2 = fnval(TWH_w2,{Pwh2,Winj2});

P1 = Pwh1;%fnval(PWH_w1,{Pwh1,Winj1});
P2 = Pwh2;%fnval(PWH_w2,{Pwh2,Winj2});

CP_G1 = fnval(FCP_G,{P1,T1});
CP_L1 = fnval(FCP_L,{P1,T1});

CP_G2 = fnval(FCP_G,{P2,T2});
CP_L2 = fnval(FCP_L,{P2,T2});

WL1 = fnval(WL_w1,{Pwh1,Winj1});
WG1 = fnval(WG_w1,{Pwh2,Winj1});


WL2 = fnval(WL_w2,{Pwh1,Winj2});
WG2 = fnval(WG_w2,{Pwh2,Winj2});

% H_Gw1 = WG1*(fnval(FH_G,{P1,T1})-fnval(FH_G,{1e5,0}));
% H_Gw2 = WG2*(fnval(FH_G,{P2,T2})-fnval(FH_G,{1e5,0}));
% 
% H_Lw1 = WL1*(fnval(FH_L,{P1,T1})-fnval(FH_L,{1e5,0}));
% H_Lw2 = WL2*(fnval(FH_L,{P2,T2})-fnval(FH_L,{1e5,0}));

H_Gw1 = CP_G1*WG1*T1;
H_Lw1 = CP_L1*WL1*T1;

H_Gw2 = CP_G2*WG2*T2;
H_Lw2 = CP_L2*WL2*T2;

Hin = H_Gw1+H_Lw1+H_Gw2+H_Lw2;

%T_INL = ((CP_G1*WG1+CP_L1*WL1)*T1+(CP_G2*WG2+CP_L2*WL2)*T2)/((CP_G1*WG1+CP_L1*WL1)+(CP_G2*WG2+CP_L2*WL2));
T_INL = Hin/((CP_G1*WG1+CP_L1*WL1)+(CP_G2*WG2+CP_L2*WL2));

alpha1 = (WG1)/(WG1+WL1);
alpha2 = (WG2)/(WG2+WL2);

alpha = (WG1+WG2)/(WG1+WG2+WL1+WL2);

Wtot = WL1 + WL2 + WG1 + WG2;

Pin = fnval(FP_INL,{T_INL,alpha,Wtot,Ptop});

%T_top = fnval(FT_TOP,{T_INL,alpha,Wtot,Ptop});
T_top = fnval(FTH_TOP,{T_INL,Hin,Ptop});

WLtot = WL1 + WL2;
end

