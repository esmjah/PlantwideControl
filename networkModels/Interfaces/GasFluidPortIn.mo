within networkModels.Interfaces;
connector GasFluidPortIn
  output SI.Pressure p(nominal = 1e5, min = 0);
  flow input SI.MassFlowRate w(nominal = 1, min = 0);
end GasFluidPortIn;
