within networkModels.networkComponents;
model manifold

  parameter manifoldParameters par;
  import SI = Modelica.SIunits;
  import Modelica.Constants.pi;
  import Modelica.Constants.R;
  import g = Modelica.Constants.g_n;

  parameter SI.Area A(fixed = true) = pi * par.r ^ 2;
  // Cross section area of pipe (m2)er (m2)
  parameter Real softPApprox(fixed = true) = 1;
  parameter Real softApprox(fixed = true) = 0.001;

  Interfaces.TwoPhaseFluidPortIn fin(p(start=2.17e6),w(start={2,29}));
  Interfaces.TwoPhaseFluidPortOut fout(p(start=5e6),w(start={-2,-29}));

//   SI.Density rho_G(fixed = false, min = 0, start = 16.75);
//   SI.Density rho_mix(fixed = false, min = 0,start = 200);
//   SI.Velocity Uslin(fixed = false, min = 0.5, max = 3, start = 1);
//   SI.Velocity Uslin_max(fixed = false, min = 0.5, max = 3, start = 1);
//   SI.Velocity Usgin(fixed = false, min = 2, max = 6, start = 4);
//   SI.Velocity Usgin_max(fixed = false, min = 2, max = 6, start = 4);
//   SI.Velocity Umin(fixed = false, min = 0.1, start = 5);
//   SI.ReynoldsNumber Re(fixed = false, min = 0, start = 1e3);
//   Real Lambda(fixed = false, min = 0, start = 0.01);
//   SI.DynamicViscosity vism(fixed = false, min = 0,start=0.00001);
  parameter SI.Pressure Fric(start = 1e5) = 18105;

equation
//   rho_G = fout.p*par.M_G/(R*par.T);
//   rho_mix = par.Alpha_L_av * par.rho_L + (1 - par.Alpha_L_av) * rho_G;
//
//   // Uslin_max = min(fin.w[2] / (A * par.rho_L) , 3)
//   // y := (-((x1-x2)^2+0.01^2)^0.5/2+(x1+x2)/2);
//   Uslin_max= (-((fin.w[2] / (A * par.rho_L)-3)^2+softApprox^2)^0.5/2+(fin.w[2] / (A * par.rho_L)+3)/2);
//   // Uslin = max(Uslin_max,0.5)
//   // y := ((x1 - x2) ^ 2 + 0.01 ^ 2) ^ 0.5 / 2 + (x1 + x2) / 2;
//   Uslin = ((Uslin_max - 0.5) ^ 2 + softApprox ^ 2) ^ 0.5 / 2 + (Uslin_max + 0.5) / 2;
//
//   //Usgin_max = min(fin.w[1]/(A * rho_G) , 6)
//   // y := (-((x1-x2)^2+0.01^2)^0.5/2+(x1+x2)/2);
//   Usgin_max = (-((fin.w[1]/(A * rho_G)-6)^2+softApprox^2)^0.5/2+(fin.w[1]/(A * rho_G)+6)/2);
//   // Usgin = max(Usgin_max,2)
//   // y := ((x1 - x2) ^ 2 + 0.01 ^ 2) ^ 0.5 / 2 + (x1 + x2) / 2;
//   Usgin = ((fin.w[1]/(A * rho_G) - 2) ^ 2 + softApprox ^ 2) ^ 0.5 / 2 + (fin.w[1]/(A * rho_G) + 2) / 2;
//
//   Umin = Uslin + Usgin;
//   vism = par.Alpha_L_av * par.visl + (1 - par.Alpha_L_av) * par.visg;
// //  Re = (rho_mix * Umin * (2 * par.r) / vism);
//   Re = sqrt((rho_mix * Umin * (2 * par.r) / vism)^2+softPApprox)/2+(rho_mix * Umin * (2 * par.r) / vism)/2;
//   Lambda=(1/(-1.8*log((par.eps/(2 * par.r)/3.7)^1.11+6.9/Re)/log(10)))^2;

  //Fric = Lambda * par.rho_L * Uslin ^ 2 * par.L / (2 * par.r);

  fout.w[1] = -fin.w[1];
  fout.w[2] = -fin.w[2];

  fin.p = fout.p + Fric;

end manifold;
