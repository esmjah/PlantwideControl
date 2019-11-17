within networkModels.networkComponents;
model PID

parameter Real Kc=1;
parameter Real Ti=600;
parameter Real Td=1;

parameter Real maxSetPoint=1;
parameter Real minSetPoint=0;

parameter Real Bias=0;
parameter Real Tf=0.1;

Real intError(start=0,nominal=1e7);

Real Error;
Real setPoint;
Real temp;
Real x(start = 0,fixed = true);

Modelica.Blocks.Interfaces.RealInput extSetPoint;
Modelica.Blocks.Interfaces.RealInput measurement;
Modelica.Blocks.Interfaces.RealOutput y(start = Bias,fixed = true);

equation
setPoint = minSetPoint + (maxSetPoint-minSetPoint) * extSetPoint;
Error = setPoint - measurement;

temp = -Kc*Error - Kc/Ti*intError - Kc*Td/Tf*Error;

der(intError) = Error;

der(x) = -1/Tf*x + Kc*Td/Tf^2 * Error;

y = temp+x+Bias;

end PID;
