within networkModels.networkComponents;
model gasliftWell

  import SI = Modelica.SIunits;
  import Modelica.Constants.pi;
  import Modelica.Constants.R;
  import g = Modelica.Constants.g_n;
  // Simplified model for riser slugging
  // Oil well added as inlet boundary condition
  // By: Esmaeil Jahanshahi
  // August 2011, NTNU, Norway
  // x1: Mass of gas in the well (m_Gw)
  // x2: Mass of liquid in the well (m_Lw)

  parameter gasliftWellParameters par;
  //parameter gasliftWellTuner tun(par=par);

  parameter Real softApprox(fixed = true) = 0.001;
  parameter Real softPApprox(fixed = true) = 1;
  parameter Real softDivision(fixed = true) = 0.01;

  Modelica.Blocks.Interfaces.RealInput GOR(min = -100, max = 100,start=0);
  Modelica.Blocks.Interfaces.RealInput z1(min = 0, max = 1,start=0.5);
  //Modelica.Blocks.Interfaces.RealInput P_res(min = 150e5, max = 200e5, start=160e5);
  Modelica.Blocks.Interfaces.RealInput d_PI(min = -10, max = 10,start=0);
  Real ORF(min = 0.5, max = 150,start=4.6689);
  Real PI(min=0.5*par.PI_nom, max=2*par.PI_nom, start=par.PI_nom);
  //SI.MassFlowRate uInjAnnulus(min = 0, max = 3,start=0.7);

//Interfaces.GasFluidPortIn inj(p(start=5.6e6,nominal=1e5),w(start=1.4,nominal=1));
  Interfaces.GasFluidPortIn fin(p(start=6.5e6, nominal=1e5),w(start=1.4,nominal=1));
  Interfaces.TwoPhaseFluidPortOut fout(p(start=25e5, nominal=1e5),w(start={-0.00001,-0.000001}));
//// **************** Annulus Parameters ************************
  SI.Pressure P_a_b(start=7.5e6);
  SI.Density rho_G_a_b(start=56);

//// ************  Calculating Friction of Riser NOMINAL **********

parameter SI.Length r_r=par.D_w/2;
parameter SI.Volume V_r=pi*((par.D_w/2)^2)*par.L_h+pi*((par.D_w/2)^2)*par.L_r;
parameter SI.Volume V_bh = par.L_bh*par.S_bh;

Real alpha_G_m_in(min=0,max=1,start=0);
//Real Z_c(min=0.5, max=1.5, start=1);

// Should this be an average?
SI.Velocity Usl_av(min=0, start=par.w_res_in/(par.rho_L*pi*(r_r^2)));
SI.Velocity Usg_av(start=10,nominal=10,min=0);
SI.Velocity U_M_av(start=11,nominal=10,min=0);

//// Density of gas inside riser
SI.Density rho_G_r(start=11,min=0);

//// ************* Average mixture density ***************
SI.Density rho_mix(min=0);
//// Pressure riser top
Real alpha_L_av(min=0,max=1);
SI.Pressure P_r_t(start=2.6e6,nominal=20e5,min=10e5);

SI.DynamicViscosity my_mix(start=1.0e-4);

SI.ReynoldsNumber Re(start=2e6,nominal=1e5,min=0);

Real lambda(min = 0,start=0.15,nominal=1e-2);
Real u(min = 1e-4,max = 1);

SI.Pressure F_riser(start=6.6e6,nominal=1e5,min=0);
SI.Pressure deltaP_valve(start=1.5e5,nominal=1e5,min=0);
//// *********** Pressure at bottom of riser ******************

//// Liquid velocity at bottom hole
SI.Velocity Ul_b(start=2.07);
//// Reynolds number at bottom-hole:
SI.ReynoldsNumber Re_b(start=5.1e5);
//// Darcy Friction Factor
Real lambda_b(min= 0,start=0.015);
//// Pressure loss due to friction from injection point to bottom-hole:
SI.Pressure F_b(start=14e5);
//// Bottom-hole pressure
SI.Pressure P_bh(start=6.13e6,min=0,max=par.P_r,nominal=1e5);
SI.Pressure P_bh_f(start=90e5,min=80e5,max=par.P_r,nominal=90e5);
/// Injection point pressure
SI.Pressure P_r_b(start=6.13e6,min=0);

//SI.Pressure P_res_f(min = 150, max = 200, start=160, nominal=160);

//// Liquid Inflow rate
SI.MassFlowRate w_r_in(start=17,min=0);
SI.MassFlowRate w_L_r_in(start=1,min=0);
SI.MassFlowRate w_G_r_in(start=0,min=0);
SI.MassFlowRate w_G_iph(start=0,min=0);
SI.MassFlowRate w_G_inj(start=1,min=0);
SI.MassFlowRate w_G_f(start=1,min=0);
//SI.MassFlowRate w_G_f2(start=1,min=0);
//// **********  Gas density in bottom of riser **************
SI.Density rho_G_r_b(start=30,min=0);

//// Alpha liquid in
Real alpha_L_in(min=0,max=1,start=0.33);

//// ******** Liquid volume fraction top of riser ****************
SI.VolumeFraction alpha_L_t_max(min=0.01,max=1,start=0.16);
SI.VolumeFraction alpha_L_t_min(min=0.01,max=1,start=0.16);
SI.VolumeFraction alpha_L_t(min=0.01,max=1,start=0.16);
SI.VolumeFraction alpha_G_t(min=0.01,max=1,start=0.8);

//// ***********  Density mixture top of riser *************************
SI.Density rho_M_t(start=300,min=0);

//// **********  Mixture flow out of riser ****************************
SI.MassFlowRate w_mix_out(start=15,min=5,max=40);

//// **********  Liquid and gas flow out of riser ****************************
SI.MassFlowRate w_G_out(start=1,min=0.5,max=3);
SI.MassFlowRate w_L_out(start=15,min=5,max=36);

//// ********* Liquid mass fraction top of riser ********************

Real alpha_L_m_t(min=0.05,max=1,start=0.8);
Real GOR_m_t(min=0.01,max=1,start=0.2);

SI.Pressure delta_p_inj(min=1e4,nominal=1e5);

  SI.Mass m_Ga(start = 4362.64,nominal = 4500, fixed = true,min=4000,max=5000);
  SI.Mass m_Gw(start = 238, nominal = 250, fixed = true, min=180,max=300);
  SI.Mass m_Lw(start = 9201.51,nominal = 9500,fixed = true, min=8000,max=11000);
//  SI.Pressure P_bth(start = 100e5,nominal = 100,fixed = true,min=0,max=200e5);

equation
//Z_c = par.Z_ca + par.Slope_z*(P_bh_f-par.P_z);
alpha_G_m_in = 0.01*(par.GOR_nom+GOR)/(1+0.01*(par.GOR_nom+GOR));
Usl_av = (1-alpha_G_m_in)*par.w_res_in/(par.rho_L*pi*(r_r^2));

//// Density of gas inside riser
//rho_G_r=m_Gw/(V_r + par.L_bh*par.S_bh - m_Lw/par.rho_L);
rho_G_r=m_Gw/(V_r - (m_Lw-par.rho_L*par.L_bh*par.S_bh)/par.rho_L);

//// ************* Average mixture density ***************
//rho_mix=(m_Gw+m_Lw-par.rho_L*par.L_bh*par.S_bh)/(V_r);
rho_mix=(m_Gw+m_Lw-par.rho_L*par.L_bh*par.S_bh)/(par.L_r*pi*r_r^2);
alpha_L_av= m_Lw/(par.rho_L*(V_bh+V_r));
//// Pressure riser top
P_r_t=par.Z_ct*rho_G_r*R*par.T_r/(par.M_G_r_t);
//P_r_t = max(fout.p,P_r_t);

//Usg_av= (alpha_G_m_in*par.w_res_in+par.w_g_a_in)/(rho_G_r*pi*(r_r^2));
Usg_av*(rho_G_r*pi*(r_r^2)) = (alpha_G_m_in*par.w_res_in+w_G_f);

U_M_av = Usl_av+Usg_av;

my_mix = alpha_L_av*par.my_L+(1-alpha_L_av)*par.my_G;

Re=(2*rho_mix*U_M_av*r_r)/my_mix;
//Re=(2*rho_mix*U_M_av*r_r)/par.my_L;

//lambda=0.0056+(0.5*(Re^(-0.32)));
//// Darcy Friction Factor

lambda=(1/(1.8*log((par.ew/par.D_w/3.7)^1.11+6.9/Re)/log(10)))^2;
//F_riser= alpha_L_av*(lambda*rho_mix*(U_M_av^2)*(par.L_r+par.L_h))/(4*r_r);
F_riser=par.K_f*lambda*rho_mix*(U_M_av^2)*(par.L_r+par.L_h)/(4*r_r);
//// *********** Pressure at bottom of riser ******************
P_r_b = P_r_t + (rho_mix*g*par.L_r) +F_riser;

//// Liquid velocity at bottom hole
Ul_b= par.w_res_in/(par.rho_L*par.S_bh);
//// Reynolds number at bottom-hole:
Re_b=(par.rho_L*Ul_b*par.D_b)/par.my_L;
//// Friction factor at bottom-hole:
//lambda_b=0.0056+(0.5*(Re_b^(-0.32)));
//// Darcy Friction Factor
lambda_b=(1/(1.8*log((par.ew/par.D_b/3.7)^1.11+6.9/Re_b)/log(10)))^2;
//// Pressure loss due to friction from injection point to bottom-hole:
F_b=(lambda_b*par.rho_L*(Ul_b^2)*(par.L_bh))/(4*r_r);
//// Bottom-hole pressure
P_bh=P_r_b+par.rho_L*g*par.L_bh +F_b;

//// Liquid Inflow rate

PI = par.PI_nom + d_PI*par.PI_nom/50;
w_r_in = PI*(par.P_r-P_bh_f);
//w_r_in=par.PI*(((((1e5*par.P_r-P_bh) - 0) ^ 2 + softPApprox) ^ 0.5 / 2 + ((1e5*par.P_r-P_bh) + 0) / 2));

w_L_r_in = (1-alpha_G_m_in)*w_r_in;
w_G_r_in = alpha_G_m_in*w_r_in;
//// **********  Gas density in bottom of riser **************
rho_G_r_b=P_r_b*par.M_G_r_t/(R*par.T_r);

// ********* Pressure at bottom of annulus *********************
fin.p = par.Z_ca * R * par.T_a * m_Ga / (par.M_G_a * par.V_a);
P_a_b = fin.p + m_Ga * g * par.L_a / par.V_a;
// ********** Injected gas into annulus ******************
//fin.w = w_G_a_in = par.K_s*u*sqrt(rho_G_in*max(par.P_gs-P_a_t,0));
// This is now an input
// ********** Density of gas in bottom of annulus **********
rho_G_a_b = P_a_b * par.M_G_a / (R * par.T_a);

//// ************ Injected gas at bottom of well **********************
//w_G_inj = tun.K_a * sqrt(max((rho_G_a_b *(P_a_b - P_r_b)),0));
delta_p_inj = P_a_b - P_r_b;
w_G_inj = par.k1*par.K_a * sqrt((((rho_G_a_b *(delta_p_inj) - 0) ^ 2 + softPApprox ^ 2) ^ 0.5 / 2 + (rho_G_a_b *(delta_p_inj) + 0) / 2));

//// Alpha liquid in
//alpha_G = m_Gw/(m_gw+m_Lw);

alpha_L_in*(w_L_r_in*rho_G_r_b+(w_G_inj+w_G_r_in)*par.rho_L)= (w_L_r_in*rho_G_r_b);
//alpha_L_in= (w_L_r_in*rho_G_r_b)/(w_L_r_in*rho_G_r_b+(inj.w+w_G_r_in)*par.rho_L+softDivision);

//// ******** Liquid volume fraction top of riser ****************
alpha_L_t_min = par.K_alpha*2*alpha_L_av - alpha_L_in;
//alpha_L_t_max = max(alpha_L_t_min,0);
alpha_L_t_max = ((alpha_L_t_min ^ 2 + softApprox^2)^0.5)/2 + alpha_L_t_min/2;
//alpha_L_t= min(alpha_L_t_max,1);
alpha_L_t = (-((alpha_L_t_max-1)^2 + softApprox^2)^0.5)/2 + (alpha_L_t_max+1)/2;
rho_M_t= alpha_L_t*par.rho_L +(1-alpha_L_t)*rho_G_r;

//// ***********  Density mixture top of riser *************************
//rho_M_t=max(0.1,rho_M_t);

//// **********  Mixture flow out of riser ****************************

// w_mix_out=tun.K_r*z1*sqrt(rho_M_t*(max(P_r_t-fout.p,0)));
// w_mix_out=tun.K_r*z1*sqrt(rho_M_t*(((((P_r_t-fout.p) - 0) ^ 2 + softPApprox) ^ 0.5 / 2 + ((P_r_t-fout.p) + 0) / 2)));

  // y := ((x1 - x2) ^ 2 + 0.01 ^ 2) ^ 0.5 / 2 + (x1 + x2) / 2;
  u = (sqrt(z1 ^ 2 + softApprox^2)) / 2 + z1/2;
  //ORF = 1/(z1^2*par.CD^2) - 1;
  ORF = (((1/(u^2*par.k2^2*par.CD^2) -1) ^ 2 + softApprox^2) ^ 0.5) / 2 + (1/(u^2*par.k2^2*par.CD^2) -1)/2;
  //ORF = sqrt((1/(z1^2*par.CD^2) -1)^2);

  w_mix_out = par.A_c*sqrt(2*rho_M_t * ((((P_r_t - fout.p) ^ 2 + softPApprox) ^ 0.5) / 2 + (P_r_t - fout.p) / 2)/ORF);
  //deltaP_valve = (ORF*w_mix_out^2)/(2*rho_M_t*par.A_c^2);
  deltaP_valve = (P_r_t - fout.p);

/// ********* Adiabatic Flash Evaporation **************
alpha_G_t = (1-alpha_L_t)*(P_bh_f/par.P_z)^(-1/par.Gamma);

//// ********* Liquid mass fraction top of riser ********************

alpha_L_m_t*((1-alpha_G_t)*par.rho_L + alpha_G_t*rho_G_r)=((1-alpha_G_t)*par.rho_L);

//// Liquid mass flow rate out of riser
w_L_out = alpha_L_m_t*w_mix_out;
fout.w[2]= -w_L_out;

//// gas mass flow rate out of riser
w_G_out = (1-alpha_L_m_t)*w_mix_out;
fout.w[1]= -w_G_out;

// inter-phase mass transport
w_G_iph = par.K_i*w_L_r_in;

GOR_m_t = w_G_out/w_L_out;

//// Derivatives
der(w_G_f) = (fin.w - w_G_f)/par.Tp1;  // Transportation lag in annulus
der(P_bh_f) = (P_bh - P_bh_f)/par.Tp2; // Transportation lag in tubing
der(m_Ga) = w_G_f  - w_G_inj;  // Annulus volume integrating state
der(m_Gw) = w_G_r_in + w_G_iph + w_G_inj + fout.w[1];
der(m_Lw) = w_L_r_in - w_G_iph + fout.w[2];
//der(P_bth) = (P_bh - P_bth)/tun.Tau;

  annotation (DymolaStoredErrors);
end gasliftWell;
