within networkModels.Interfaces;
connector ThreePhaseFluidPortOut
  input SI.Pressure p(nominal = 1e5, min = 0);
  flow output SI.MassFlowRate w[3](each nominal = 10, each max = 0);
end ThreePhaseFluidPortOut;
