within networkModels.networkComponents;
record gasliftWellParameters

  import SI = Modelica.SIunits;
  import Modelica.Constants.pi;
  import g = Modelica.Constants.g_n;

parameter SI.Temperature T_a(min = 315, max = 350,nominal = 100) = 273+50; //348;
parameter SI.MolarMass M_G_a(min = 0.016, max = 0.025,nominal = 0.020) = 0.023; //0.01631;
parameter SI.Length L_a = 2048;
parameter SI.Length D_a = 0.2;
  parameter SI.Volume V_a = 2098*pi*(0.1^2) + 10*pi*(0.05^2); //par.L_a*pi*((par.D_a/2)^2);
parameter SI.Pressure P_gs = 140e5;
//parameter Real k1 = 1.1;
//parameter Real K_a = 0.000155590989907951;
parameter Real K_a = 1.5559e-04;
parameter Real K_f = 0.72;
parameter Real K_i = 0.00067;//0.0002;
parameter Real Gamma = 1.7576;//Adiabatic Flash Evaporation
parameter Real k1(min=0.7,max=1.5)= 1.2;   // Annulus
parameter Real k2=0.9;  // Tubing
parameter Real K_alpha(min=0.5, max=1.5) = 1;

parameter SI.Time Tp1 = 600; // Transportation lag 1
parameter SI.Time Tp2 = 900;  // Transportation lag 2

parameter Real Slope_z(fixed=true) = 0.0014e-5;//0.0014;
parameter SI.Pressure P_z(start=100e5,min=80e5,max=155e5,nominal=1e5) =  100e5;
parameter Real Z_ca = 1;
parameter Real Z_ct = 1;
// Annulus
parameter SI.MassFlowRate w_g_a_in=1.012;
//parameter Real z2 = 1.5;
parameter SI.Pressure Pat = 86.3e5;

//% Viscosity (Liquid)
parameter SI.DynamicViscosity my_L = 0.364e-3; //Read from olga values at 100 C
parameter SI.DynamicViscosity my_G = 1.5e-5;

//% Density Liquid [kg/m^3]
parameter SI.Density rho_L(min = 740, max = 800,nominal = 100)=770;//760; //Read from Olga at values at 100 C

//% Temperature riser
//par.T_r=373; //[K] avg temperature in the middle of the tubing olga
parameter SI.Temperature T_r(min = 340, max = 380,nominal = 100)=273+96;//375; //[K]temperature at top of the tubing olga

//% Gas molar weight in top of riser        [kg/mol]
parameter SI.MolarMass M_G_r_t(min = 0.016, max = 0.025,nominal = 0.020) = 0.02; //0.0167;
//% Gas molar weight in top of riser        [kg/mol]
parameter SI.MolarMass M_G_r(min = 0.016, max = 0.025,nominal = 0.020) = 0.02; //0.0167;

parameter Real PI_nom = 2.47e-6;  //Productivity Index

//% OLGA TUBING VOLUMES CALCULATION
parameter SI.Length D_w = 0.124;

parameter SI.Length L_r=2048;//vertical
parameter SI.Length L_h = 25; //Horizontal

//% Roughness of well
parameter SI.Length ew = 4.5e-5;
parameter SI.Length ep = 3e-5;

//% OLGA BOTTOM HOLE VOLUMES CALCULATION
parameter SI.Length L_bh=100;
parameter SI.Area S_bh=pi*((0.2/2)^2);

parameter SI.Length D_b = 0.2;

//% Reservoir properties
parameter SI.Pressure P_r = 160e5;
parameter Real GOR_nom(min=0, max=100) = 0;

//% Steady State Values From OLGA
parameter SI.MassFlowRate w_res_in = 17;

parameter Real z1=0.7;
parameter SI.Pressure Pbh=8728310;
//parameter SI.Pressure Pab=103.165e5;
parameter SI.Pressure Ptt=1850900;
parameter SI.Pressure Ptb=8161975;

//parameter SI.Pressure Pat=92.2582e5;

parameter Real CD = 0.84;
//parameter SI.Length D_v=0.025;
//parameter SI.Area A_v=pi*(D_v/2)^2;

parameter SI.Length D_c = 0.07;
parameter SI.Area A_c = pi*(0.07/2)^2;

//parameter Real k3=0.96897;//0.88578;  // Injection to annulus
//parameter Real K_s = k3*w_g_a_in/(z2*sqrt(rho_G_in*(140*1e5-Pat)));     //Into annulus,

end gasliftWellParameters;
