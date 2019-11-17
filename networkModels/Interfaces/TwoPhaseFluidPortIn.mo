within networkModels.Interfaces;
connector TwoPhaseFluidPortIn
  output SI.Pressure p(nominal = 1e5, min = 0);
  flow input SI.MassFlowRate w[2](each nominal = 1, each min = 0);
end TwoPhaseFluidPortIn;
