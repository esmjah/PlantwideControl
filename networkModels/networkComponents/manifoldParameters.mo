within networkModels.networkComponents;
record manifoldParameters

  import SI = Modelica.SIunits;
  import Modelica.Constants.pi;
  //////////////
  // ********** Constants and parameters for New simple model ***********
  parameter SI.Density rho_L(min = 700, max = 830) = 735;
  // Liquid density (Kg/m3)
  // TODO check units
  // Feed pipe inclination (Rad) = 1 Deg
  parameter SI.Radius r(min = 0.05, max = 0.2) = 8*0.0254/2;
  // Radius of pipe (m)
  parameter SI.Length L(min = 100, max = 400) = 800;
  // Length of pipe from well to manifold (m) - used as a tuning
  parameter SI.Temperature T(min = 320, max = 370,nominal = 100) = 273+95;
  // Pipeline temprature (K)
  parameter SI.MolarMass M_G(min = 0.020, max = 0.025) = 0.019; //0.02143;
  // Molecular weight of Gas (kg/kmol)
  //Boundary Condition in ss
  // *************************************************************************
  // TODO: what is this?
  parameter SI.DynamicViscosity visl(min = 1.2260e-4, max = 1.6260e-4) = 1.4260e-4;

  parameter SI.DynamicViscosity visg = 1.39e-5;

  parameter SI.Length eps(min = 2.6e-5, max = 3.0e-5) = 2.8e-5;
  // Roughness of pipe

  parameter Real Alpha_L_av(min = 0.2, max = 0.4) = 0.2576; // 0.257937;

end manifoldParameters;
