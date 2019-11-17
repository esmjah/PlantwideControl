within networkModels.networkComponents;
model gasManifold
  import SI = Modelica.SIunits;
  Interfaces.GasFluidPortOut fout(w(start=-1,min=-3,max=-0.5),p(start=6.8e6));
  //parameter SI.MassFlowRate constantOutput = 1;
  Modelica.Blocks.Interfaces.RealInput u;
equation
  fout.w = -u;
end gasManifold;
