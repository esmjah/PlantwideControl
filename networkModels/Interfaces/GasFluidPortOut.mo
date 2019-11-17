within networkModels.Interfaces;
connector GasFluidPortOut
  input SI.Pressure p(nominal = 1e5, min = 0);
  flow output SI.MassFlowRate w(nominal = 1, max = 0);
end GasFluidPortOut;
