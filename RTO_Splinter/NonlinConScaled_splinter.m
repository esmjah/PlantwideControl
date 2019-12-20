function C = NonlinConScaled_splinter( x,w1,w2,FP_INL,FT_TOP,FCP_G,FCP_L)
%NonlinCon Summary of this function goes here
%   Detailed explanation goes here

global P_sep DeltaP1_manifold DeltaP2_manifold temp_offset
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

Pwh1 = x(6) + Delta_Pwh1 + DeltaP1_manifold;
Pwh2 = x(6) + Delta_Pwh2 + DeltaP2_manifold;
Ptop = P_sep + Delta_Ptop;


T1 = TWH_w1.eval([Pwh1 Winj1]);
T2 = TWH_w2.eval([Pwh2 Winj2]);

P1 = x(6);%PWH_w1.eval([Pwh1 Winj1]);
P2 = x(6);%PWH_w2.eval([Pwh2 Winj2]);

CP_G1 = FCP_G.eval([P1 T1]);
CP_L1 = FCP_L.eval([P1 T1]);

CP_G2 = FCP_G.eval([P2 T2]);
CP_L2 = FCP_L.eval([P2 T2]);

WL1 = WL_w1.eval([Pwh1 Winj1]);
WG1 = WG_w1.eval([Pwh1 Winj1]);


WL2 = WL_w2.eval([Pwh2 Winj2]);
WG2 = WG_w2.eval([Pwh2 Winj2]);

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
T_INL = Hin/((CP_G1*WG1+CP_L1*WL1)+(CP_G2*WG2+CP_L2*WL2)) - temp_offset;

alpha = 100*(WG1+WG2)/(WG1+WG2+WL1+WL2);

Wtot = WL1 + WL2 + WG1 + WG2;

%Hin = Saturate(Hin,3.6361e+06,8.4804e+06);
Ptop = Saturate(Ptop,5.1,6.2);
T_INL = Saturate(T_INL,93,97);
Wtot = Saturate(Wtot,25,37);
alpha = Saturate(alpha,5,9);

Pm = FP_INL.eval([T_INL alpha Wtot Ptop]);

T_top = FT_TOP.eval([T_INL alpha Wtot Ptop]);
%T_top = FTH_TOP.eval([T_INL Hin Ptop]);

ceq = x(6) - Pm;
c = T_top - 10;

C = [ceq;c];

end

