within networkModels.Interfaces;
connector TwoPhaseFluidPortOut
  input SI.Pressure p(nominal = 1e5, min = 0);
  flow output SI.MassFlowRate w[2](each nominal = 10, each max = 0);
end TwoPhaseFluidPortOut;
