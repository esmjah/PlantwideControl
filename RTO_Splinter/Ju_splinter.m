function Jx = Ju_splinter( x,w1,w2,FP_INL,FCP_G,FCP_L,Lambda)
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

PriceL = Lambda(1);
PriceG = Lambda(2);

T1 = TWH_w1.eval([Pwh1 Winj1]);
T2 = TWH_w2.eval([Pwh2 Winj2]);

TWH_w1_jacs = TWH_w1.eval_jacobian([Pwh1 Winj1]);
dT1_du1 = TWH_w1_jacs(2);
dT1_du3 = TWH_w1_jacs(1);

TWH_w2_jacs = TWH_w2.eval_jacobian([Pwh2 Winj2]);
dT2_du2 = TWH_w2_jacs(2);
dT2_du4 = TWH_w1_jacs(1);

P1 = Pwh1;%PWH_w1.eval([Pwh1 Winj1]);
P2 = Pwh2;%PWH_w2.eval([Pwh2 Winj2]);

PWH_w1_jacs = PWH_w1.eval_jacobian([Pwh1 Winj1]);
dP1_du1 = PWH_w1_jacs(2);

PWH_w2_jacs = PWH_w2.eval_jacobian([Pwh2 Winj2]);
dP2_du2 = PWH_w2_jacs(2);


CP_G1 = FCP_G.eval([P1 T1]);
CP_L1 = FCP_L.eval([P1 T1]);

CP_G1_jacs = FCP_G.eval_jacobian([P1 T1]);
CP_L1_jacs = FCP_L.eval_jacobian([P1 T1]);

dCP_G1_dP1 = CP_G1_jacs(1);
dCP_G1_dT1 = CP_G1_jacs(2);

dCP_L1_dP1 = CP_L1_jacs(1);
dCP_L1_dT1 = CP_L1_jacs(2);

CP_G2 = FCP_G.eval([P2 T2]);
CP_L2 = FCP_L.eval([P2 T2]);

CP_G2_jacs = FCP_G.eval_jacobian([P2 T2]);
CP_L2_jacs = FCP_L.eval_jacobian([P2 T2]);

dCP_G2_dP2 = CP_G2_jacs(1);
dCP_G2_dT2 = CP_G2_jacs(2);

dCP_L2_dP2 = CP_L2_jacs(1);
dCP_L2_dT2 = CP_L2_jacs(2);

WL1 = WL_w1.eval([Pwh1 Winj1]);
WG1 = WG_w1.eval([Pwh1 Winj1]);

WL1_jacs = WL_w1.eval_jacobian([Pwh1 Winj1]);
WG1_jacs = WG_w1.eval_jacobian([Pwh1 Winj1]);

dWL1_du1 = WL1_jacs(2);
dWL1_du3 = WL1_jacs(1);

dWG1_du1 = WG1_jacs(2);
dWG1_du3 = WG1_jacs(1);

WL2 = WL_w2.eval([Pwh2 Winj2]);
WG2 = WG_w2.eval([Pwh2 Winj2]);

WL2_jacs = WL_w2.eval_jacobian([Pwh2 Winj2]);
WG2_jacs = WG_w2.eval_jacobian([Pwh2 Winj2]);

dWL2_du2 = WL2_jacs(2);
dWL2_du4 = WL2_jacs(1);

dWG2_du2 = WG2_jacs(2);
dWG2_du4 = WG2_jacs(1);

% H_Gw1 = WG1*(fnval(FH_G,{P1,T1})-fnval(FH_G,{1e5,0}));
% H_Gw2 = WG2*(fnval(FH_G,{P2,T2})-fnval(FH_G,{1e5,0}));
% 
% H_Lw1 = WL1*(fnval(FH_L,{P1,T1})-fnval(FH_L,{1e5,0}));
% H_Lw2 = WL2*(fnval(FH_L,{P2,T2})-fnval(FH_L,{1e5,0}));

H_Gw1 = CP_G1*WG1*T1;
H_Lw1 = CP_L1*WL1*T1;

dCP_G1_du1 = dCP_G1_dP1*dP1_du1 + dCP_G1_dT1*dT1_du1;
dCP_L1_du1 = dCP_L1_dP1*dP1_du1 + dCP_L1_dT1*dT1_du1;

dH_Gw1_du1 = dCP_G1_du1*(WG1*T1) + dWG1_du1*(CP_G1*T1) + dT1_du1*(CP_G1*WG1);
dH_Lw1_du1 = dCP_L1_du1*(WL1*T1) + dWL1_du1*(CP_L1*T1) + dT1_du1*(CP_L1*WL1);

H_Gw2 = CP_G2*WG2*T2;
H_Lw2 = CP_L2*WL2*T2;

dCP_G2_du2 = dCP_G2_dP2*dP2_du2 + dCP_G2_dT2*dT2_du2;
dCP_L2_du2 = dCP_L2_dP2*dP2_du2 + dCP_L2_dT2*dT2_du2;

dH_Gw2_du2 = dCP_G2_du2*(WG2*T2) + dWG2_du2*(CP_G2*T2) + dT2_du2*(CP_G2*WG2);
dH_Lw2_du2 = dCP_L2_du2*(WL2*T2) + dWL2_du2*(CP_L2*T2) + dT2_du2*(CP_L2*WL1);

Hin = H_Gw1+H_Lw1+H_Gw2+H_Lw2;

dHin_du1 = dH_Gw1_du1 + dH_Lw1_du1;
dHin_du2 = dH_Gw2_du2 + dH_Lw2_du2;

T_DEN = (CP_G1*WG1+CP_L1*WL1)+(CP_G2*WG2+CP_L2*WL2);

dT_DEN_du1 = dCP_G1_du1*WG1 + dWG1_du1*CP_G1 + dCP_L1_du1*WL1 + dWL1_du1*CP_L1;
dT_DEN_du2 = dCP_G2_du2*WG2 + dWG2_du2*CP_G2 + dCP_L2_du2*WL2 + dWL2_du2*CP_L2;

%T_INL = ((CP_G1*WG1+CP_L1*WL1)*T1+(CP_G2*WG2+CP_L2*WL2)*T2)/((CP_G1*WG1+CP_L1*WL1)+(CP_G2*WG2+CP_L2*WL2));
T_INL = Hin/T_DEN - temp_offset;

dT_INL_du1 = (dHin_du1*T_DEN - dT_DEN_du1*Hin)/(T_DEN^2);
dT_INL_du2 = (dHin_du2*T_DEN - dT_DEN_du2*Hin)/(T_DEN^2);

alpha = 100*(WG1+WG2)/(WG1+WG2+WL1+WL2);

Wtot = WL1 + WL2 + WG1 + WG2;
dWtot_du1 = dWL1_du1 + dWG1_du1;
dWtot_du2 = dWL2_du2 + dWG2_du2;

Ptop = Saturate(Ptop,5.1,6.2);
T_INL = Saturate(T_INL,93,97);
Wtot = Saturate(Wtot,22,37);
alpha = Saturate(alpha,5,9);

%disp(FP_INL.eval_jacobian([T_INL alpha Wtot Ptop]))
P_INL_jacs = 0*FP_INL.eval_jacobian([T_INL alpha Wtot Ptop]);
dP_INL_dT_INL = P_INL_jacs(1);
dP_INL_dAlpha = P_INL_jacs(2);
dP_INL_dWtot = P_INL_jacs(3);
%dP_INL_dPtot = P_INL_jacs(4);

dAlpha_du1 = 100*(dWG1_du1*Wtot - dWtot_du1*(WG1+WG2))/(Wtot^2);
dAlpha_du2 = 100*(dWG2_du2*Wtot - dWtot_du2*(WG1+WG2))/(Wtot^2);

dPm_du1 = dP_INL_dT_INL*dT_INL_du1 + dP_INL_dAlpha*dAlpha_du1 + dP_INL_dWtot*dWtot_du1;
dPm_du2 = dP_INL_dT_INL*dT_INL_du2 + dP_INL_dAlpha*dAlpha_du2 + dP_INL_dWtot*dWtot_du2;

Jx1 = PriceG - PriceL*(dWL1_du1 + dWL1_du3*dPm_du1 + dWL2_du4*dPm_du1);
Jx2 = PriceG - PriceL*(dWL2_du2 + dWL2_du4*dPm_du2 + dWL1_du3*dPm_du2);
Jx3 = - PriceL*dWL1_du3;
Jx4 = - PriceL*dWL2_du4;
Jx5 = - PriceL*dWL1_du3 - PriceL*dWL2_du4;
Jx6 = - PriceL*dWL1_du3 - PriceL*dWL2_du4;

Jx = [Jx1;Jx2;Jx3;Jx4;Jx5;Jx6];

end

