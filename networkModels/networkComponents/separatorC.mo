within networkModels.networkComponents;
model separatorC
  import SI = Modelica.SIunits;
  Interfaces.ThreePhaseFluidPortIn fin(p(start=constantPressure));
  parameter SI.Pressure constantPressure = 50e5;
equation
  fin.p = constantPressure;
end separatorC;
