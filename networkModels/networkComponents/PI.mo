within networkModels.networkComponents;
model PI

parameter Real Kc=1;
parameter Real Ti=600;

parameter Real maxSetPoint=1;
parameter Real minSetPoint=0;

parameter Real softApprox = 1e-5;
parameter Real satMax = 1e10;
parameter Real satMin = -1e10;

Real intError(start=0,nominal=1e7);

Real error;
Real setPoint;

Modelica.Blocks.Interfaces.RealInput extSetPoint;
Modelica.Blocks.Interfaces.RealInput measurement;
Modelica.Blocks.Interfaces.RealOutput u(start = 0,fixed = true);
Real u0;
Real uSatMax;

equation
setPoint = minSetPoint + (maxSetPoint-minSetPoint) * extSetPoint;
error = setPoint - measurement;

u0 = -Kc*error - Kc/Ti*intError;

uSatMax = ((u0 - satMin) ^ 2 + softApprox) ^ 0.5 / 2 + (u0 + satMin) / 2;

u = (-((uSatMax-satMax)^2+softApprox)^0.5/2+(uSatMax+satMax)/2);

der(intError) = error;

end PI;
