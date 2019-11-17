within networkModels.networkComponents;
model separatorTwoPhase
  import SI = Modelica.SIunits;
  Interfaces.TwoPhaseFluidPortIn fin;
  //parameter SI.Pressure constantPressure = 50e5;
  parameter SI.Pressure Psep_nom = 510000;
  Modelica.Blocks.Interfaces.RealInput d_Psep(min=-2e5,max=2e5,start=0);
  SI.Pressure Psep(min=4e5,max=6e5,start=Psep_nom);
equation
  Psep = Psep_nom + d_Psep;
  fin.p = Psep;
end separatorTwoPhase;
