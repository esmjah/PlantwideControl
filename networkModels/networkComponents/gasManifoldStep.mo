within networkModels.networkComponents;
model gasManifoldStep

  import SI = Modelica.SIunits;
  Interfaces.GasFluidPortOut fout(w(start=-1,min=-3,max=-0.5),p(start=6.8e6));
  parameter SI.MassFlowRate offset = 1;
  parameter SI.MassFlowRate height = 0;
  parameter SI.Time startTime = 0;
  SI.MassFlowRate flowRate(start=1,min=0.5,max=3);
equation
  flowRate = offset + (if time < startTime then 0 else height);
  fout.w = -flowRate;

end gasManifoldStep;
