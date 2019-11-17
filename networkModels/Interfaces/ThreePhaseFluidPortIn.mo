within networkModels.Interfaces;
connector ThreePhaseFluidPortIn
  output SI.Pressure p(nominal = 1e5, min = 0);
  flow input SI.MassFlowRate w[3](each nominal = 1, each min = 0);
end ThreePhaseFluidPortIn;
