within networkModels.networkComponents;
model wellheadTwoPhaseVarLiquid
  import SI = Modelica.SIunits;
  Interfaces.TwoPhaseFluidPortOut fout;
  Modelica.Blocks.Interfaces.RealInput massFlowIn;
  parameter Real glrIn = 0.0416;
  Real w_g_in = massFlowIn * glrIn;
  Real w_l_in = massFlowIn * (1 - glrIn);
  /// TODO: correct, left as Thomas did
equation
  fout.w = -{w_g_in, w_l_in};
end wellheadTwoPhaseVarLiquid;
