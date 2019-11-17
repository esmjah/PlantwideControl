within networkModels.networkComponents;
model twoPhaseSource
  import SI = Modelica.SIunits;
  Interfaces.TwoPhaseFluidPortOut fout(p(start=25e6),w(start={-wG,-wL}));
  parameter SI.MassFlowRate wG = 1;
  parameter SI.MassFlowRate wL = 30;
equation
  fout.w[1] = -wG;
  fout.w[2] = -wL;
end twoPhaseSource;
