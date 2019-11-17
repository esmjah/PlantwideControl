within networkModels.networkComponents;
model wellheadTwoPhase
  import SI = Modelica.SIunits;
  Interfaces.TwoPhaseFluidPortOut fout;
  parameter SI.MassFlowRate constantOutput[2] = {1, 1};
equation
  fout.w = -constantOutput;
end wellheadTwoPhase;
