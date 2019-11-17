within networkModels.networkComponents;
record pipelineRiserParemeters
  import SI = Modelica.SIunits;
  import Modelica.Constants.pi;
  import Modelica.Constants.R;
  import g = Modelica.Constants.g_n;
  //////////////
   // *******************************************************
  // *************** Tuning Parameters *********************
   parameter Real k_h(min = 0.5, max = 1.2, nominal=1) = 0.74;//0.65;
   parameter Real K_gT(min = 0.5, max = 2, nominal=1) = 0.83; //0.9
   parameter Real K_oT(min = 0.5, max = 4, nominal=1) = 2; //2.33
   parameter Real K_alpha(min = 0.5, max = 1.5) = 0.72;//0.73
   parameter Real K_i = 0.014;
   parameter Real Gamma = 1.7576;//Adiabatic Flash Evaporation
   parameter Real K_f = 0.7;//0.4

  // **********************************************************
  parameter Real Slope_z(fixed=true) = -0.00143;
  parameter SI.Pressure P_z = 2.10134e+006;
  parameter Real Z_c0 = 0.95;
  // **************** Physical Parameters *******************
  parameter SI.Density rho_L(min = 700, max = 830) = 750;//733.4;
  // Liquid density (Kg/m3)
  parameter SI.Density rho_w(min = 900, max = 1100) = 1000;
  // TODO check units
  parameter SI.Angle theta(min = 0.5 * pi / 180, max = 2 * pi / 180) = pi / 180;
  // Feed pipe inclination (Rad) = 1 Deg
  parameter SI.Radius r1(min = 0.05, max = 0.2) = 8*0.0254/2;
  // Radius of pipe (m)
  parameter SI.Radius r2(min = 0.04, max = 0.2) = 8*0.0254/2;
  // Radius of riser (m)
  parameter SI.Length Lm(min = 100, max = 400) = 700;
  // Length of pipe from well to manifold (m) - used as a tuning
  parameter SI.Length L1(min = 4100, max = 4400) = 4300;
  // Length of upstream pipe (m)
  parameter SI.Height L2(min = 280, max = 320) = 300;
  // Height of riser
  parameter SI.Length L3(min = 90, max = 110) = 100;
  // length of horozontal top section (m)
  parameter SI.Temperature T1(min = 320, max = 370,nominal = 100) = 273+95;
  // Pipeline temprature (K)
  parameter SI.Temperature T2(min = 320, max = 370,nominal = 100) = 273+71.5;// 298.3;
  // Riser temprature (K)
  // TODO: discuss units of M_G with Esmaeil
  parameter SI.MolarMass M_Gp(min = 0.020, max = 0.025) = 0.020; //0.02143;
  parameter SI.MolarMass M_Gr(min = 0.020, max = 0.025) = 0.020; //0.02143;
  parameter SI.Time Tf = 400;  // Transportation lag

  // Molecular weight of Gas (kg/kmol)
  //Boundary Condition in ss
  // *************************************************************************
  parameter SI.DynamicViscosity visl(min = 1.2260e-4, max = 1.6260e-4) = 1.4260e-4;

  parameter SI.DynamicViscosity visg = 1.39e-5;

  parameter SI.Length eps(min = 2.6e-5, max = 3.0e-5) = 2.8e-5;
  // Roughness of pipe

  parameter Real Cd = 0.83; // Choke valve discharge coefficient
  parameter SI.Length D_c = 7*0.0254;  // Choke valve dimatere
  parameter SI.Area A_c = pi*(7*0.0254/2)^2;  // Choke valve dimatere
  //parameter Real Alpha_L1_ss(min = 0.2, max = 0.4) = 0.2576; // 0.257937;
  parameter SI.VolumeFraction Alpha_L1_ss = 0.2576; //wL_in * rho_G1 / (wL_in * rho_G1 + wG_in * rho_L); 0.2576;
  // ****************************************************************
  // ******** Steady State Values Baesd on OLGA *************
   parameter Real z0(min = 0.03, max = 1) = 0.5;
   parameter SI.Pressure P0(min = 2.0e5, max = 70.2e5) = 5.156e5;
   // Pressure after choke valve (Pa)
   parameter SI.Pressure P1(min = 10e5, max = 40e5) = 22.7905e5;
   parameter SI.Pressure P2_t(min = 5e5, max = 7e5) = 5.96184e5;
   parameter SI.MassFlowRate wG_in(min = 0, max = 100) = 2.104;
   // TODO: check units!
   parameter SI.MassFlowRate wL_in(min = 0, max = 100) = 29.8761;
  //********************************************************************
  //*********************** calculated parameters  ********************
  //parameter Real K_gT = 0.0223328;
  //parameter Real K_oT = 0.174109;

  parameter Real K_g = v_G1 / sqrt((P1 - Fric_pipe - P2_b) / rho_G1);
  parameter Real K_o = wL_in / (A_l * sqrt(rho_L * (P1 - Fric_pipe + rho_L * g * h1ss - P2_b)));

  //parameter Real K_pc = K_pcT * w_mix / (z0 * sqrt(rho_t * (P2_t - P0)));
protected
  parameter Real softApprox(fixed = true) = 0.001;
  parameter SI.Area A1 = pi * r1 ^ 2;
  // Cross section area of pipeline (m2)
  parameter SI.Area A2 = pi * r2 ^ 2;
  // Cross section area of riser (m2)
  parameter SI.Height hc = 2 * r1 / cos(theta);
  // Critical liquid level (m)
  parameter SI.MassFlowRate w_mix = wL_in + wG_in;
  parameter SI.Density rho_G1 = P1 * M_Gp / (R * T1);
  parameter SI.Height h1ss = k_h * Alpha_L1_ss * hc;
  parameter SI.Density rho_mix1 = Alpha_L1_ss*rho_L + (1-Alpha_L1_ss)*rho_G1;
  parameter SI.Velocity Uslin = wL_in /(A1 * rho_L);
  parameter SI.Velocity Usgin = wG_in/(A1*rho_G1);
  parameter SI.Velocity Umin = Uslin + Usgin;
  parameter SI.DynamicViscosity vism_p = Alpha_L1_ss*visl + (1-Alpha_L1_ss)*visg;
  parameter SI.ReynoldsNumber Re1 = rho_mix1*Umin*(2*r1)/vism_p;
  parameter Real Lambda1 = 0.0056 + 0.5 * Re1 ^ (-0.32);

  parameter SI.Pressure Fric_pipe = 0.5 * Alpha_L1_ss * Lambda1 * rho_L * Uslin ^ 2 * L1 / (2 * r1);
  parameter SI.Pressure P2_av = (P1 + rho_L * g * h1ss - Fric_pipe + P2_t) / 2;
  parameter SI.Density rho_G2_av = P2_av * M_Gr / (R * T2);
  parameter SI.MassFraction Alpha_L2_av = rho_G2_av * wL_in / (rho_G2_av * wL_in + rho_L * wG_in);
  parameter SI.Density rho_mix_av = rho_G2_av * (1 - Alpha_L2_av) + rho_L * Alpha_L2_av;
  parameter SI.Velocity Usl2 = wL_in / (rho_L * A2);
  parameter SI.Velocity Usg2 = wG_in / (rho_G2_av * A2);
  parameter SI.Velocity Um = Usl2 + Usg2;
  parameter SI.DynamicViscosity vism_r = Alpha_L2_av*visl + (1-Alpha_L2_av)*visg;
  parameter SI.ReynoldsNumber Re2 = rho_mix_av * Um * 2 * r2 / vism_r;

  parameter Real Lambda2 = (1/(-1.8*log((eps/(2*r2)/3.7)^1.11+6.9/Re2)/log(10)))^2;
  parameter SI.Pressure Fric_riser = 0.5*Lambda2*rho_mix_av*Um^2*(L2+L3)/(2*r2);
  parameter SI.Pressure P2_b = P2_t + rho_mix_av * g * L2 + Fric_riser;
  parameter SI.Density rho_G2 = P2_t * M_Gr / (R * T2);
  parameter SI.MassFraction Alpha_Lt = rho_G2 * wL_in / (rho_G2 * wL_in + rho_L * wG_in);
  parameter SI.Density rho_t = Alpha_Lt * rho_L + (1 - Alpha_Lt) * rho_G2;
  //parameter SI.Area A_g = A1 * (maxFunc(((hc - h1ss)/hc),0))^2;
  parameter SI.Area A_g = A1 * ((((hc - h1ss) / hc - 0) ^ 2 + softApprox ^ 2) ^ 0.5 / 2 + ((hc - h1ss) / hc + 0) / 2) ^ 2;
  parameter SI.Area A_l = A1 - A_g;
  parameter SI.Velocity v_G1 = wG_in / (rho_G1 * A_g);

end pipelineRiserParemeters;
