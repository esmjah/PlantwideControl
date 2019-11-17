within networkModels.networkComponents;
model firstOrder

  import SI = Modelica.SIunits;
  parameter Real k(unit="1")=1 "Gain";
  parameter SI.Time T(start=1) "Time Constant";

  Modelica.Blocks.Interfaces.RealInput u;
  Modelica.Blocks.Interfaces.RealOutput y;

//initial equation
//    der(y) = 0;
equation

  der(y) = (k*u - y)/T;

end firstOrder;
