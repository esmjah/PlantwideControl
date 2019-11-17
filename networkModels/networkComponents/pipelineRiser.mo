within networkModels.networkComponents;
model pipelineRiser
  // Simplified model for riser slugging
  // By: Esmaeil Jahanshahi
  // August 2009, NTNU, Norway
  parameter networkModels.networkComponents.pipelineRiserParemeters par;
  import SI = Modelica.SIunits;
  import Modelica.Constants.pi;
  import Modelica.Constants.R;
  import g = Modelica.Constants.g_n;
  //  parameters derivated from constants and dimensions
  //
  //  input: flow & pressure
  //  calcualate -> densities alphas velocities reynolds lambdas
  /// output -> DP
  // geometry
  parameter SI.Area A1(fixed = true) = pi * par.r1 ^ 2;
  // Cross section area of pipeline (m2)
  parameter SI.Area A2(fixed = true) = pi * par.r2 ^ 2;
  // Cross section area of riser (m2)
  parameter SI.Volume V2(fixed = true) = A2 * (par.L2 + par.L3);
  // Total volume in riser (m3)
  parameter SI.Volume V1(fixed = true) = A1 * par.L1;
  // Total volume in pipeline (m3)
  parameter SI.Height hc(fixed = true) = 2 * par.r1 / cos(par.theta);
  // Critical liquid level (m)
  parameter Real softApprox(fixed = true) = 0.001;
  parameter Real softPApprox(fixed = true) = 1;
  Interfaces.TwoPhaseFluidPortIn fin(p(start=2.17e6,nominal=1e5),w(start={2,29}));
  Interfaces.TwoPhaseFluidPortOut fout(p(start=6e5,nominal=1e5),w(start={-2,-29}));
  Modelica.Blocks.Interfaces.RealInput z(min = 0, max = 1,start=0.5);
  //Modelica.Blocks.Interfaces.RealInput Pin_ss(min = 19e5, max = 30e5,start=23e5);

  Real ORF(min = 0.5, max = 150,start=4.6689);
  SI.Height h1(fixed = false, min = 0, start = 0.02);
  SI.Density rho_G1(fixed = false, min = 0, start = 16.75);
//  SI.Density rho_ss(fixed = false, min = 0, start = 16.75);
  SI.Velocity Uslin(fixed = false, min = 0, start = 1);
  SI.ReynoldsNumber Re1(fixed = false, min = 0, start = 1e3);
  Real Lambda1(fixed = false, min = 0, start = 0.01);
  SI.Pressure Fric_pipe(fixed = false, start = 1e5);
  SI.Volume V_G2(fixed = false, min = 0, start = 11);
  SI.Pressure P2_t(fixed = false, start = 5.8e5, nominal = 1e5);
  SI.Pressure P1(fixed = false, start = par.P1, min=20e5, max=30e5, nominal = 1e5);
  SI.Pressure P1_f(fixed = false, start = par.P1, min=20e5, max=30e5, nominal = 1e5);
  SI.Density rho_G2(fixed = false, min = 0, start = 5);
  Real Alpha_L2_av(fixed = false, min = 0, max = 1, start = 0.12);
  SI.Density rho_mix_av(fixed = false, min = 0, start = 70);
  SI.Velocity Usl2(fixed = false, min = 0, start = 1);
  SI.Velocity Usg2(fixed = false, min = 0, start = 14);
  SI.Velocity Um(fixed = false, min = 0, start = 15);
  SI.ReynoldsNumber Re2(fixed = false, min = 0, start = 1e6);
  Real u(min = 1e-4,max = 1);
  Real Lambda2(fixed = false, min = 0, start = 0.01);
  SI.Pressure Fric_riser(fixed = false, min = 0, start = 0.3e5);
  SI.Area A_g(fixed = false, min = 0.002,max=A1, start = 0.03);
  SI.Area A_l(fixed = false, min = 0.002,max=A1, start = 0.006);
  SI.Pressure P2_b(fixed = false, start = 10e5);
  SI.MassFlowRate w_G1(fixed = false, min = 0, start = 2);
  SI.MassFlowRate w_L1(fixed = false, min = 0, start = 30);
  SI.MassFlowRate w_G_out(start=2,min=0.5,max=4);
  SI.MassFlowRate w_L_out(start=30,min=5,max=40);
  SI.MassFlowRate w_G_iph(start=0.02,min=0,max=1);
  SI.MassFraction Alpha_Lb(fixed = false, min = 0.05, max = 0.8, start = 0.2);
  SI.MassFraction Alpha_LbMax;
  SI.VolumeFraction Alpha_L1_av(fixed = false, min = 0.01, max = 0.8, start = 0.2576);
  SI.VolumeFraction Alpha_Lt_max(fixed = false, min = 0.01, max = 0.8, start = 0.05);
  SI.VolumeFraction Alpha_Lt_min(fixed = false, min = 0.01, max = 0.8, start = 0.05);
  SI.VolumeFraction Alpha_Lt(fixed = false, min = 0.01, max = 0.8, start = 0.05);
  SI.VolumeFraction Alpha_Gt(fixed = false, min = 0.01, max = 0.8, start = 0.05);
  SI.MassFraction Alpha_Lmt(fixed = false, min = 0.1, max = 1, start = 0.89);
  SI.Density rho_t(fixed = false, min = 0, start = 65);
  SI.MassFlowRate w_mix_out(fixed = false, min = 3, max=40, start = 29);
  SI.Density rho_mix1(fixed = false, min = 0,start = 200);
  SI.DynamicViscosity vism_p(fixed = false, min = 0,start=0.00001);
  SI.Velocity Usgin(fixed = false, min = 0,start=0.00001);
  SI.Velocity Umin(fixed = false, min = 0,start=0.00001);
  SI.DynamicViscosity vism_r(fixed = false, min = 1e-7,start=0.00001,nominal=1e-5);
  SI.Mass m_lp0t;
  SI.Height h1ss;
  Real Z_c(min=0.5, max=1.5, start=1);

  SI.Pressure deltaP_G(min = 0.1e5,nominal=1e5);
  SI.Pressure deltaP_L(min = 0.1e5,nominal=1e5);
  SI.Pressure deltaP_valve(start=0.5e5,nominal=1e5,min=0);

  // Iniital Conditions
  //******************//***********************//
  //               States & Variables                         //
  //******************//***********************//

  //  parameter SI.Mass m_lp0t(fixed=false,start=25229.3714534857);

  /// liquid level in steady-state
  SI.Mass m_gp(start = 1466.35, fixed = true, min = 1100,max=1600,nominal= 1280);
  // mass of gas in the horizontal pipeline
  SI.Mass m_lp(start = 26353.6, fixed = true, min = 26000,max=26460,nominal= 26300);
  // mass of oil in the horizontal pipeline
  SI.Mass m_gr(start = 43.8007, fixed = true, min = 35, max = 50,nominal= 48);
  // mass of gas in the riser
  SI.Mass m_lr(start = 2724.82, fixed = true, min = 1700,max=3500,nominal= 1900);

equation
  fin.p = P1_f;
  Z_c = par.Z_c0; // + 1e-5*(P1_f-par.P_z)*par.Slope_z; //P1_f/par.P1;;
  rho_G1 = fin.p*par.M_Gp/(Z_c*R*par.T1); //m_gp / (V1 - m_lp / par.rho_L);
  Alpha_L1_av = fin.w[2] * rho_G1 / (fin.w[2] * rho_G1 + fin.w[1] * par.rho_L);
  m_lp0t = V1 * par.rho_L * par.Alpha_L1_ss;
  //// some tunning here
  //  parameter SI.Height h1ss(fixed=false,min=0,start=0.0553652169975689);     /// liquid level in steady-state
  h1ss = par.k_h * Alpha_L1_av * hc;
  //h1 = max(h1ss + sin(par.theta) / (A1 * (1 - Alpha_L1_av) * par.rho_L) * (m_lp - m_lp0t),0);
  h1 = ((h1ss + sin(par.theta) / (A1 * (1 - Alpha_L1_av) * par.rho_L) * (m_lp - m_lp0t) - 0) ^ 2 + softApprox) ^ 0.5 / 2 + (h1ss + sin(par.theta) / (A1 * (1 - Alpha_L1_av) * par.rho_L) * (m_lp - m_lp0t) + 0) / 2;
  // h1 = min(h1,40)
  //h1 = (-((h0-0.2)^2+0.01^2)^0.5/2+(h0+0.2)/2);
  //h1 = h1ss;

  V_G2 = V2 - m_lr / par.rho_L;
  P1 = Z_c * m_gp * R * par.T1 / (par.M_Gp * (V1 - m_lp / par.rho_L));

  rho_mix1 = Alpha_L1_av * par.rho_L + (1 - Alpha_L1_av) * rho_G1;
  Uslin = fin.w[2] / (A1 * par.rho_L);
  Usgin*(A1 * rho_G1) = fin.w[1];
  Umin = Uslin + Usgin;
  vism_p = Alpha_L1_av * par.visl + (1 - Alpha_L1_av) * par.visg;
//  Re1 = (rho_mix1 * Umin * (2 * par.r1) / vism_p);
  Re1 = sqrt((rho_mix1 * Umin * (2 * par.r1) / vism_p)^2+softPApprox)/2+(rho_mix1 * Umin * (2 * par.r1) / vism_p)/2;
  //Lambda1=(1/(-1.8*log((par.eps/(2 * par.r1)/3.7)^1.11+6.9/Re1)/log(10)))^2;//  WHAT LAMBDA TO USE?

  //Lambda1 = 0.0056 + 0.5 * Re1 ^ (-0.32);
  Lambda1 = 0.0056 + 0.5 * (((Re1 - softApprox) ^ 2 + softApprox ^ 2) ^ 0.5 / 2 + (Re1 + softApprox) / 2) ^ (-0.32);
  Fric_pipe = Alpha_L1_av * Lambda1 * par.rho_L * Uslin ^ 2 * par.L1 / (2 * par.r1);

  //P1 =  P1_f - Lambda1 * par.rho_L * Uslin ^ 2 * par.Lm / (2 * par.r1);
  //Fric_pipe = Alpha_L1_av * Lambda1 * par.rho_L * Uslin ^ 2 * par.L1 / (2 * par.r1);

  //P2_t = max(m_gr*R*par.T2/(par.M_Gr*V_G2),fout.p);

  P2_t = sqrt((Z_c*m_gr*R*par.T2/(par.M_Gr*V_G2) - fout.p) ^ 2 + softPApprox ^ 2) / 2 + (Z_c*m_gr*R*par.T2/(par.M_Gr*V_G2) + fout.p) / 2;
  //P2_t*(par.M_Gr * V_G2) = ((m_gr * R * par.T2 - (par.M_Gr * V_G2)*fout.p) ^ 2 + softPApprox) ^ 0.5 / 2 + (m_gr * R * par.T2+ (par.M_Gr * V_G2)*fout.p) / 2;
  rho_G2*V_G2 = m_gr;
  Alpha_L2_av = m_lr / (par.rho_L * V2);
  rho_mix_av = (m_gr + m_lr) / V2;
  Usl2 = fin.w[2] / (par.rho_L * A2);
  Usg2*(rho_G2 * A2) = fin.w[1];
  Um = Usl2 + Usg2;
  vism_r = Alpha_L2_av * par.visl + (1 - Alpha_L2_av) * par.visg;

  //Re2 = (rho_mix_av * Um * (2 * par.r2)/vism_r);
  Re2 = sqrt((rho_mix_av * Um * (2 * par.r2)/vism_r)^2+softPApprox)/2+(rho_mix_av * Um * (2 * par.r2)/vism_r)/2;
  //Lambda2=0.0056 + 0.5*Re2^(-0.32);
  Lambda2 = 0.0056 + 0.5 * (((Re2 - softApprox) ^ 2 + softApprox ^ 2) ^ 0.5 / 2 + (Re2 + softApprox) / 2) ^ (-0.32);
 // Lambda2 = (1 / (-1.8 * log((par.eps / (2*par.r2) / 3.7) ^ 1.11 + 6.9 / Re2) / log(10))) ^ 2;
 // Lambda2 = (1 / (-1.8 * log((par.eps / (2*par.r2) / 3.7) ^ 1.11 + 6.9 / ((Re2^2 + softApprox^2)^0.5/2+Re2/2)) / log(10))) ^ 2;
  Fric_riser = Alpha_L2_av * Lambda2 * rho_mix_av * Um ^ 2 * (par.L2 + par.L3) / (2 * par.r2);
  //Fric_riser = par.K_f * Lambda2 * rho_mix_av * Um ^ 2 * (par.L2 + par.L3) / (2 * par.r2);
  //A_g = (h1<par.hc)*(A1/par.hc^2)*(par.hc - h1)^2;
  //A_g = A1 * (maxFunc(((hc - h1)/hc),0))^2;
  A_g = A1 * ((((hc - h1) / hc - 0) ^ 2 + softApprox) ^ 0.5 / 2 + ((hc - h1) / hc + 0) / 2) ^ 2;
  A_l = A1 - A_g;
  P2_b = P2_t + rho_mix_av * g * par.L2 + Fric_riser;
  //w_G1 = par.K_gT*par.K_g*A_g*sqrt(rho_G1*max  (0,P1_f-Fric_pipe-P2_b));
  //w_G1 = par.K_gT*par.K_g*A_g*sqrt(rho_G1*maxFunc(P1_f-Fric_pipe-P2_b,0));
  deltaP_G = P1_f - Fric_pipe - P2_b;
  w_G1 = par.K_gT * par.K_g * A_g * sqrt(rho_G1 * (((deltaP_G) ^ 2 + softApprox) ^ 0.5 / 2 + (deltaP_G) / 2));
  //w_L1 = (P1_f-Fric_pipe+par.rho_L*g*h1-P2_b>0)*...
  //       par.K_oT*par.K_o*A_l*sqrt(par.rho_L*abs(P1_f-Fric_pipe+par.rho_L*g*h1-P2_b));
  //w_L1 = par.K_oT * par.K_o*A_l*sqrt(par.rho_L*(P1_f-Fric_pipe+par.rho_L*g*h1-P2_b));
  deltaP_L = P1_f - Fric_pipe + par.rho_L * g * h1 - P2_b;
  w_L1 = par.K_oT * par.K_o * A_l * sqrt(par.rho_L * (((deltaP_L) ^ 2 + softApprox) ^ 0.5 / 2 + (deltaP_L) / 2));

  //Alpha_Lb = 1 -  A_g/A1;
  //if(Alpha_Lb<=Alpha_L2_av)
  //    Alpha_Lb=Alpha_L2_av;
  //end
  //Alpha_Lb = maxFunc(1 -  A_g/A1,Alpha_L2_av);
  Alpha_LbMax = ((1 - A_g / A1 - Alpha_L2_av) ^ 2 + softApprox) ^ 0.5 / 2 + (1 - A_g / A1 + Alpha_L2_av) / 2;
  //Alpha_Lb = minFunc(Alpha_LbMax,1);
  Alpha_Lb =(-((Alpha_LbMax-1)^2+softApprox)^0.5/2+(Alpha_LbMax+1)/2);
  //Alpha_Lt =  2*Alpha_L2_av - Alpha_Lb;
  //if(Alpha_Lt>Alpha_L2_av)
  //    Alpha_Lt = Alpha_L2_av;
  //elseif(Alpha_Lt<0)
  //    Alpha_Lt = 0;
  //end
  //Alpha_Lt =  maxFunc(minFunc(2*Alpha_L2_av - Alpha_Lb,Alpha_L2_av),0);
  Alpha_Lt_min = par.K_alpha*2*Alpha_L2_av - Alpha_Lb;
  //alpha_Lt_max = max(alpha_Lt_min,0);
  Alpha_Lt_max = ((Alpha_Lt_min ^ 2 + softApprox^2)^0.5)/2 + Alpha_Lt_min/2;
  //alpha_Lt= min(alpha_Lt_max,Alpha_L2_av);
  Alpha_Lt = (-((Alpha_Lt_max-Alpha_L2_av)^2 + softApprox^2)^0.5)/2 + (Alpha_Lt_max+Alpha_L2_av)/2;

  /// ********* Adiabatic Flash Evaporation **************
  Alpha_Gt = (1-Alpha_Lt); //*(P1_f/par.P_z)^(-1/par.Gamma);

  //Alpha_Lt = (((-((2 * Alpha_L2_av - Alpha_Lb - Alpha_L2_av) ^ 2 + softApprox) ^ 0.5 / 2) + (2 * Alpha_L2_av - Alpha_Lb + Alpha_L2_av) / 2 - 0) ^ 2 + softApprox) ^ 0.5 / 2 + ((-((2 * Alpha_L2_av - Alpha_Lb - Alpha_L2_av) ^ 2 + softApprox) ^ 0.5 / 2) + (2 * Alpha_L2_av - Alpha_Lb + Alpha_L2_av) / 2 + 0) / 2;
  Alpha_Lmt*((1-Alpha_Gt) * par.rho_L + Alpha_Gt * rho_G2) = (1-Alpha_Gt) * par.rho_L;
  rho_t = Alpha_Lt * par.rho_L + (1 - Alpha_Lt) * rho_G2;
  // if(z>0.95 && par.Cd = 1)
  //     par.Cd = 0.95;
  // end

  // y := ((x1 - x2) ^ 2 + 0.01 ^ 2) ^ 0.5 / 2 + (x1 + x2) / 2;
  u = (sqrt(z ^ 2 + softApprox^2)) / 2 + z/2;
  //ORF = 1/(z^2*par.Cd^2) - 1;
  ORF = (sqrt((1/(u^2*par.Cd^2) -1) ^ 2 + softApprox^2)) / 2 + (1/(u^2*par.Cd^2) -1)/2;
  //ORF = sqrt((1/(z^2*par.Cd^2) -1) ^ 2);
  w_mix_out = par.A_c* sqrt(2*rho_t * ((((P2_t - fout.p) ^ 2 + softPApprox)^0.5)/2 + (P2_t - fout.p)/2)/ORF);
  //deltaP_valve = (ORF*w_mix_out^2)/(2*rho_t*par.A_c^2);
  deltaP_valve = (P2_t - fout.p);
  //w_mix_out = A2*sqrt(2*rho_t*max(0,P2_t-fout.p)/ORF);
  //w_mix_out = par.K_pc*z*sqrt(rho_t*max(0,P2_t-fout.p));
  //w_mix_out = par.K_pc*z*sqrt(rho_t*maxFunc(0,(P2_t-fout.p)));
  //w_mix_out = par.K_pc * z * sqrt(rho_t * (((0 - (P2_t - fout.p)) ^ 2 + softApprox) ^ 0.5 / 2 + (0 + P2_t - fout.p) / 2));

  // inter-phase mass transport
  w_G_iph = par.K_i*w_L1;

  w_G_out = (1-Alpha_Lmt)*w_mix_out;
  w_L_out = w_mix_out - w_G_out;

  fout.w[1] = - w_G_out;
  fout.w[2] = - w_L_out;

  // State equattions
  der(P1_f) = (P1 - P1_f)/par.Tf;
  der(m_gp) = fin.w[1] - w_G1;
  der(m_lp) = fin.w[2] - w_L1;
  der(m_gr) = w_G1 + w_G_iph + fout.w[1];
  der(m_lr) = w_L1 - w_G_iph + fout.w[2];
  //der(Alpha_Lin) =  (rho_ss*fin.w[2]/(rho_ss*fin.w[2]+par.rho_L*fin.w[1]) - Alpha_Lin)/par.Tau_alpha;
  annotation (DymolaStoredErrors);
end pipelineRiser;
