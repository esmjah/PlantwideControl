within networkModels.networkComponents;
model CSTR

  import SI = Modelica.SIunits;
  import Modelica.Constants.pi;
  import Modelica.Constants.R;
  import g = Modelica.Constants.g_n;

  parameter SI.Time Tau(min = 40, max = 80,nominal = 60) = 60;
  parameter SI.MolecularConcentration C_A0(min=0, max=1, nominal=1, start=0.498) = 0.498;
  parameter SI.MolecularConcentration C_B0(min=0, max=1, nominal=1, start=0.502) = 0.502;
  parameter SI.Temperature T0(min=0, max=300, nominal=600, start=426.803) = 426.803;

  Modelica.Blocks.Interfaces.RealInput T_i(min = 300, max = 600,start=424.292);
  Modelica.Blocks.Interfaces.RealInput C_Ai(min = 0, max = 1,start=1);
  Modelica.Blocks.Interfaces.RealInput C_Bi(min = 0, max = 1,start=0);

  SI.MolecularConcentration C_A(min=0, max=1, nominal=1, start=0.498);
  SI.MolecularConcentration C_B(min=0, max=1, nominal=1, start=0.502);
  SI.Temperature T(min=0, max=300, nominal=600, start=426.803);

  Real r(min=0,max=1);

equation
  r = 5000*exp(-10000/(1987*T))*C_A - 10^6*exp(-15000/(1987*T))*C_B;

  der(C_A) = (C_Ai*C_A0/0.498 - C_A)/Tau - r;
  der(C_B) = (C_Bi*C_B0/0.502 - C_B)/Tau + r;
  der(T)   = (T_i*T0/426.803 - T)/Tau  + 5*r;

  annotation (DymolaStoredErrors);
end CSTR;
