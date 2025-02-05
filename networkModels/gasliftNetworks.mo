within networkModels;
package gasliftNetworks

  model well1
    replaceable package components = networkModels.networkComponents;
    import SI = Modelica.SIunits;
    import Modelica.Constants.R;

    parameter Real PCA1InitialOutput = 0.5;
    parameter Real[3] x0 = { 3576.58, 250.926, 7702.3};

    networkModels.networkComponents.gasliftWell_Pres
                           w1(par=networkModels.Instances.well1(),    m_Ga(start=x0[1]),m_Gw(start=x0[2]),m_Lw(start=x0[3]));

    parameter SI.Pressure PCA1InitialSetPoint =  R * w1.par.T_a * x0[1] / (w1.par.M_G_a * w1.V_a);
    components.separatorTwoPhase sep(constantPressure = 2.3234e+006,fin(w(start={2.7,24.3})));

    components.PI PCA1(Kc=1e-5,Ti=600,Bias=PCA1InitialOutput,maxSetPoint=110e5,minSetPoint=70e5,satMax=1,satMin=0);
    components.PI PCT1(Kc=-0.004e-5,Ti=900,Bias=(PCA1InitialSetPoint - PCA1.minSetPoint)  / (PCA1.maxSetPoint - PCA1.minSetPoint),satMax=1,satMin=0);

    components.firstOrder PCT1Filter(k=1,T=900,y(start=(PCA1InitialSetPoint - PCA1.minSetPoint)  / (PCA1.maxSetPoint - PCA1.minSetPoint)));

    Modelica.Blocks.Interfaces.RealInput whd1PressureSetPoint;

  equation
    connect(w1.fout, sep.fin);
    connect(whd1PressureSetPoint,PCT1.extSetPoint);
    PCT1.measurement = w1.P_r_t;

    connect(PCT1.u,PCT1Filter.u);
    connect(PCT1Filter.y,PCA1.extSetPoint);
    PCA1.measurement = w1.fin.p;

    connect(PCA1.u,w1.z1);

  end well1;

  partial model well2
    replaceable package components = networkModels.networkComponents;
      networkModels.networkComponents.gasliftWell_Pres
                             w2(par=networkModels.Instances.well2(),    m_Ga(start=3896),m_Gw(start=386.42),m_Lw(start=4247.52));
  equation

  end well2;

  model wellSim "This will simulate a well"
    replaceable package components = networkModels.networkComponents;

    well1 pidWell;

    Modelica.Blocks.Sources.Constant u1(k = 2.55069e+006);

    components.gasManifold gm(constantOutput = 1,fout(p(start=86.5e5)));

  equation
    connect(gm.fout, pidWell.w1.fin);
    connect(u1.y, pidWell.whd1PressureSetPoint);

    annotation (experiment(
        StopTime=72000,
        Interval=1,
        Tolerance=1e-012), __Dymola_experimentSetupOutput);
  end wellSim;

  partial model netPres

    import SI = Modelica.SIunits;

    replaceable package components = networkModels.networkComponents;

    parameter Real[17] x0 = {160, 1.29636, 9.99259e+006, 4854.88, 337.866, 8070.42, 170, 1.32454, 1.04266e+007, 5027.63, 330.052, 8553.26, 2.5434e+006, 1812.43, 26909, 46.4096, 1563.75};

    // x0 = {160, 1, 1.01741e+007, 4697.47, 277.545, 8937.33, 170, 1, 1.06498e+007, 4845.08, 270.267, 9505.04, 2.18215e+006, 1554.94, 26912, 45.5072, 1802.1};

    SI.Pressure Dp_w1(min=0,nominal=1e5);
    SI.Pressure Dp_w2(min=0,nominal=1e5);

    networkModels.networkComponents.gasliftWell_Pres
                           w1(par=networkModels.Instances.well1(P_z=x0[3]),
                              P_res(start=x0[1], nominal=x0[1]),
                              w_G_f(start=x0[2], nominal=x0[2]),
                              P_bh_f(start=x0[3], nominal=x0[3]),
                              m_Ga(start=x0[4], nominal=x0[4]),
                              m_Gw(start=x0[5], nominal=x0[5]),
                              m_Lw(start=x0[6], nominal=x0[6]));

    networkModels.networkComponents.gasliftWell_Pres
                           w2(par=networkModels.Instances.well2(P_z=x0[9]),
                              P_res(start=x0[7], nominal=x0[7]),
                              w_G_f(start=x0[8], nominal=x0[8]),
                              P_bh_f(start=x0[9], nominal=x0[9]),
                              m_Ga(start=x0[10], nominal=x0[10]),
                              m_Gw(start=x0[11], nominal=x0[11]),
                              m_Lw(start=x0[12], nominal=x0[12]));

    components.manifold m( Fric=18105);

    components.pipelineRiser p(P1_f(start=x0[13], nominal=x0[13]),
                               m_gp(start=x0[14], nominal=x0[14]),
                               m_lp(start=x0[15], nominal=x0[15]),
                               m_gr(start=x0[16], nominal=x0[16]),
                               m_lr(start=x0[17], nominal=x0[17]));

    components.separatorTwoPhase sep(Psep(start=sep.Psep_nom),fin(w(start={2.7,36.3})));

  equation
    Dp_w1 = w1.P_r_t - p.fin.p;
    Dp_w2 = w2.P_r_t - p.fin.p;

    connect(w1.fout, m.fin);
    connect(w2.fout, m.fin);
    connect(m.fout, p.fin);
    connect(p.fout, sep.fin);
  end netPres;

  model netOpenLoopSim
    extends networkModels.gasliftNetworks.netPres;
    Modelica.Blocks.Sources.Constant u1(k = 0.5);
    components.gasManifold gm1(u = 1,fout(p(start=86.5e5)));
    Modelica.Blocks.Sources.Constant u2(k = 0.5);
    components.gasManifold gm2(u = 1,fout(p(start=91.6e5)));
    Modelica.Blocks.Sources.Constant up(k = 0.5);

    Modelica.Blocks.Sources.Constant Pr1(k = 160);
    Modelica.Blocks.Sources.Constant Pr2(k = 170);

    Modelica.Blocks.Sources.Constant GOR1(k = 0);
    Modelica.Blocks.Sources.Constant GOR2(k = 0);

    Modelica.Blocks.Sources.Constant d_PI1(k = 0);
    Modelica.Blocks.Sources.Constant d_PI2(k = 0);

    Modelica.Blocks.Sources.Constant d_Psep(k = 0);

  equation
    connect(GOR1.y, w1.GOR);
    connect(GOR2.y, w2.GOR);
    connect(d_PI1.y, w1.d_PI);
    connect(d_PI2.y, w2.d_PI);
    connect(Pr1.y, w1.P_res);
    connect(Pr2.y, w2.P_res);
    connect(gm1.fout, w1.fin);
    connect(gm2.fout, w2.fin);
    connect(u1.y, w1.z1);
    connect(u2.y, w2.z1);
    connect(up.y,p.z);
    connect(d_Psep.y,sep.d_Psep)
    annotation (experiment(
        StopTime=36000,
        Interval=1,
        Tolerance=1e-006,
        Algorithm="Dassl"),
        __Dymola_experimentSetupOutput);

  end netOpenLoopSim;

  model pidNet
    extends networkModels.gasliftNetworks.netPres;

    import SI = Modelica.SIunits;
    import Modelica.Constants.R;

  parameter Real u0PCA1(min = 0, max = 1) = 0.500782;
  parameter Real y0PCA1(min = 60e5, max = 100e5) = R * w1.par.T_a * x0[3] / (w1.par.M_G_a * w1.par.V_a);
  parameter Real KcPCA1 = 1e-5;
  parameter Real TiPCA1 = 600;
  parameter Real x0PCA1 = -TiPCA1*u0PCA1/KcPCA1;

  components.PI PCA1(Kc=KcPCA1,Ti=TiPCA1, maxSetPoint=110e5, minSetPoint=70e5, satMax=1, satMin=0, intError(start=x0PCA1, nominal=3.00469e7), u(start=u0PCA1));

  parameter Real u0PCT1(min = 0, max = 1) = (y0PCA1 - PCA1.minSetPoint)  / (PCA1.maxSetPoint - PCA1.minSetPoint);
  parameter Real y0PCT1(min = 2e6, max = 3e6) = x0[4]*R*w1.par.T_r/((w1.V_r - (x0[5]-w1.par.rho_L*w1.par.L_bh*w1.par.S_bh)/w1.par.rho_L)*w1.par.M_G_r_t);
  Modelica.Blocks.Interfaces.RealInput whd1PressureSetPoint(start=y0PCT1);
  parameter Real KcPCT1 = -0.005e-5;
  parameter Real TiPCT1 = 900;
  parameter Real x0PCT1 = -TiPCT1*u0PCT1/KcPCT1;

  components.PI PCT1(Kc=KcPCT1,Ti=TiPCT1, satMax=1,satMin=0,intError(start=x0PCT1, nominal=9.511e9), u(start = u0PCT1));

  parameter Real u0PCA2(min = 0, max = 1) = 0.500306;
  parameter Real y0PCA2(min = 60e5, max = 100e5) = R * w2.par.T_a * x0[8] / (w2.par.M_G_a * w2.par.V_a);
  parameter Real KcPCA2 = 1e-5;
  parameter Real TiPCA2 = 600;
  parameter Real x0PCA2 = -TiPCA2*u0PCA2/KcPCA2;

  components.PI PCA2(Kc=KcPCA2,Ti=TiPCA2, maxSetPoint=115e5, minSetPoint=75e5, satMax=1, satMin=0, intError(start=x0PCA2, nominal=3.00184e7), u(start=u0PCA2));

  parameter Real u0PCT2(min = 0, max = 1) = (y0PCA2 - PCA2.minSetPoint)  / (PCA2.maxSetPoint - PCA2.minSetPoint);
  parameter Real y0PCT2(min = 2e6, max = 3e6) = x0[9]*R*w2.par.T_r/((w2.V_r - (x0[10]-w2.par.rho_L*w2.par.L_bh*w2.par.S_bh)/w2.par.rho_L)*w2.par.M_G_r_t);
  Modelica.Blocks.Interfaces.RealInput whd2PressureSetPoint(start=y0PCT2);
  parameter Real KcPCT2 = -0.005e-5;
  parameter Real TiPCT2 = 900;
  parameter Real x0PCT2 = -TiPCT2*u0PCT2/KcPCT2;

  components.PI PCT2(Kc=-0.004e-5,Ti=900, satMax=1, satMin=0, intError(start=x0PCT2, nominal=9.04405e9), u(start = u0PCT2));

  parameter Real u0PCINL(min = 0, max = 1) = 0.500059;
  parameter Real y0PCINL(min = 2e6, max = 3e6) = x0[11] * R * p.par.T1 / (p.par.M_Gp * (p.V1 - x0[12] / p.par.rho_L));
  Modelica.Blocks.Interfaces.RealInput pipelinePressureSetPoint(start=y0PCINL);
  parameter Real KcPCINL = 0.2e-5;
  parameter Real TiPCINL = 900;
  parameter Real x0PCINL = -TiPCINL*u0PCINL/KcPCINL;

  components.PI PCINL(Kc=KcPCINL, Ti=TiPCINL, satMax=1,satMin=0,intError(start=x0PCINL, nominal=2.25027e8), u(start = u0PCINL));

  components.firstOrder PCT1Filter(k=1, T=600, y(start = u0PCT1, nominal = 0.422711));
  components.firstOrder PCT2Filter(k=1, T=600, y(start = u0PCT2, nominal = 0.401958));

  equation
  connect(whd1PressureSetPoint,PCT1.extSetPoint);
  PCT1.measurement = w1.P_r_t;

  connect(whd2PressureSetPoint,PCT2.extSetPoint);
  PCT2.measurement = w2.P_r_t;

  connect(PCT1.u,PCT1Filter.u);
  connect(PCT1Filter.y,PCA1.extSetPoint);
  PCA1.measurement = w1.fin.p;

  connect(PCT2.u,PCT2Filter.u);
  connect(PCT2Filter.y,PCA2.extSetPoint);
  PCA2.measurement = w2.fin.p;

  connect(PCA1.u,w1.z1);
  connect(PCA2.u,w2.z1);

  connect(pipelinePressureSetPoint,PCINL.extSetPoint);
  PCINL.measurement = p.fin.p;

  connect(PCINL.u,p.z);
  //connect(pipelinePressureSetPoint,p.Pin_ss);

    annotation (experiment(
        StopTime=36000,
        Interval=10,
        Tolerance=1e-010,
        Algorithm="Dassl"), __Dymola_experimentSetupOutput);
  end pidNet;

  model pidNetSim

    replaceable package components = networkModels.networkComponents;

    import SI = Modelica.SIunits;
    import Modelica.Constants.R;

    pidNet netPid;

    //parameter Real[3] PresssureSetpoints = { 2.3244e+006, 2.55069e+006, 2.56983e+006};

    //2.54381000e+06;
    parameter SI.Pressure PCT1SetPoint= netPid.x0[4]*R*netPid.w1.par.T_r/((netPid.w1.V_r - (netPid.x0[5]-netPid.w1.par.rho_L*netPid.w1.par.L_bh*netPid.w1.par.S_bh)/netPid.w1.par.rho_L)*netPid.w1.par.M_G_r_t);

    //2.57117000e+06;
    parameter SI.Pressure PCT2SetPoint= netPid.x0[9]*R*netPid.w2.par.T_r/((netPid.w2.V_r - (netPid.x0[10]-netPid.w2.par.rho_L*netPid.w2.par.L_bh*netPid.w2.par.S_bh)/netPid.w2.par.rho_L)*netPid.w2.par.M_G_r_t);

    //2.27915000e+06;
    parameter SI.Pressure PCINLSetPoint = netPid.x0[11] * R * netPid.p.par.T1 / (netPid.p.par.M_Gp * (netPid.p.V1 - netPid.x0[12] / netPid.p.par.rho_L));

     components.gasManifold gm1(u = 1,fout(p(start=86.5e5)));
     components.gasManifold gm2(u = 1,fout(p(start=91.6e5)));

     Modelica.Blocks.Sources.Constant up(k = PCINLSetPoint);

     Modelica.Blocks.Sources.Constant u1(k = PCT1SetPoint);
     Modelica.Blocks.Sources.Constant u2(k = PCT2SetPoint);

  equation
     connect(gm1.fout, netPid.w1.fin);
     connect(gm2.fout, netPid.w2.fin);

     connect(u1.y, netPid.whd1PressureSetPoint);
     connect(u2.y, netPid.whd2PressureSetPoint);
     connect(up.y, netPid.pipelinePressureSetPoint);

    annotation (experiment(
        StopTime=72000,
        Interval=1,
        Tolerance=1e-012,
        Algorithm="Dassl"), __Dymola_experimentSetupOutput);
  end pidNetSim;

  model pipeSim

    replaceable package components = networkModels.networkComponents;
    import SI = Modelica.SIunits;
    import Modelica.Constants.R;

    parameter Real[4] x0 = {1731.54, 27071.6, 52.4886, 839.07};

    Modelica.Blocks.Sources.Constant up(k = 0.5);

    networkModels.networkComponents.twoPhaseSource  src(wG=2.10413,wL=29.8761,fout(p(start=22.7e5)));

    networkModels.networkComponents.pipelineRiser
                                     p(
                                     m_gp(start=x0[1]),
                                     m_lp(start=x0[2]),
                                     m_gr(start=x0[3]),
                                     m_lr(start=x0[4]));

    components.separatorTwoPhase sep(constantPressure = 515600,fin(w(start={20,36.3})));

    parameter Real P1ss = 22.7905e5;// x0[1] * R * p.par.T1 / (p.par.M_G * (p.V1 - x0[2] / p.par.rho_L)); //
    Modelica.Blocks.Sources.Constant pipelineSteadyPressure(k = P1ss);

  equation
    connect(src.fout, p.fin);
    connect(p.fout, sep.fin);
    connect(up.y,p.z);
    connect(pipelineSteadyPressure.y,p.Pin_ss);
    annotation (experiment(
        StopTime=72000,
        Interval=1,
        Tolerance=1e-006), __Dymola_experimentSetupOutput);
  end pipeSim;

  model netPres_PIDstruc2
    extends networkModels.gasliftNetworks.netPres;

    import SI = Modelica.SIunits;
    import Modelica.Constants.R;
    import Modelica.Constants.pi;

    parameter Real u0PCINL(min = 0, max = 1) = 0.573572;
    parameter Real ZcINL = p.par.Z_c0; //x0[13]/p.par.P1; //1 + 1e-5*(x0[13] - p.par.P_z)*p.par.Slope_z;
    parameter Real y0PCINL(min = 2e6, max = 4e6) = ZcINL * x0[14] * R * p.par.T1 / (p.par.M_Gp * (p.V1 - x0[15] / p.par.rho_L));
    //Modelica.Blocks.Interfaces.RealInput pipelinePressureSetPoint(start=y0PCINL);
    parameter Real KcPCINL = 2e-6;
    parameter Real TiPCINL = 300;
    parameter Real x0PCINL = -TiPCINL*u0PCINL/KcPCINL;

    components.PI PCINL(Kc=KcPCINL, Ti=TiPCINL, maxSetPoint=40e5, minSetPoint=20e5, satMax=1,satMin=0,intError(start=x0PCINL, nominal=3.00778e7), u(start = u0PCINL));

    parameter Real u0PCTOP(min = 0, max = 1) = (y0PCINL - PCINL.minSetPoint)  / (PCINL.maxSetPoint - PCINL.minSetPoint);
    parameter Real ZcTOP = p.par.Z_c0; //x0[13]/p.par.P1; //1 + 1e-5*(x0[13] - p.par.P_z)*p.par.Slope_z;
    parameter Real y0PCTOP(min = 0) = ZcTOP*x0[16]*R*p.par.T2/(p.par.M_Gr*((pi * p.par.r2 ^ 2) * (p.par.L2 + p.par.L3) - x0[17] / p.par.rho_L)) - sep.Psep_nom;
    Modelica.Blocks.Interfaces.RealInput valveDP_top_SetPoint(start=y0PCTOP);
    parameter Real KcPCTOP = -2e-7;
    parameter Real TiPCTOP = 300;
    parameter Real x0PCTOP = -TiPCTOP*u0PCTOP/KcPCTOP;

    components.PI PCTOP(Kc=KcPCTOP, Ti=TiPCTOP, satMax=1,satMin=0,intError(start=x0PCTOP, nominal=1.08894e8), u(start = u0PCTOP));

    parameter Real u0PCA1(min = 0, max = 1) = 0.580034;
    parameter Real ZcA1 = w1.par.Z_ca; //x0[3]/100e5;//1 + 1e-5*(x0[3] - w1.par.P_z)*w1.par.Slope_z;
    parameter Real y0PCA1(min = 60e5, max = 100e5) = ZcA1 * R * w1.par.T_a * x0[4] / (w1.par.M_G_a * w1.par.V_a);
    parameter Real KcPCA1 = 3e-5;
    parameter Real TiPCA1 = 600;
    parameter Real x0PCA1 = -TiPCA1*u0PCA1/KcPCA1;

    components.PI PCA1(Kc=KcPCA1,Ti=TiPCA1, maxSetPoint=110e5, minSetPoint=70e5, satMax=1, satMin=0, intError(start=x0PCA1, nominal=3.02824e7), u(start=u0PCA1));

    parameter Real u0PCT1(min = 0, max = 1) = (y0PCA1 - PCA1.minSetPoint)  / (PCA1.maxSetPoint - PCA1.minSetPoint);
    parameter Real ZcT1 = w1.par.Z_ct; //x0[3]/100e5; //1 + 1e-5*(x0[3] - w1.par.P_z)*w1.par.Slope_z;
    parameter Real y0PCT1(min = 0) = ZcT1*x0[5]*R*w1.par.T_r/((w1.V_r - (x0[6]-w1.par.rho_L*w1.par.L_bh*w1.par.S_bh)/w1.par.rho_L)*w1.par.M_G_r_t) - y0PCINL - m.Fric;
    Modelica.Blocks.Interfaces.RealInput valveDP_whd1_SetPoint(start=y0PCT1);
    parameter Real KcPCT1 = -4e-9;
    parameter Real TiPCT1 = 100;
    parameter Real x0PCT1 = -TiPCT1*u0PCT1/KcPCT1;

    components.PI PCT1(Kc=KcPCT1,Ti=TiPCT1, satMax=1,satMin=0,intError(start=x0PCT1, nominal=9.25742e9), u(start = u0PCT1));

    parameter Real u0PCA2(min = 0, max = 1) = 0.607185;
    parameter Real ZcA2 = w2.par.Z_ca; //x0[9]/100e5; //1 + 1e-5*(x0[9] - w1.par.P_z)*w2.par.Slope_z;
    parameter Real y0PCA2(min = 60e5, max = 100e5) = ZcA2 * R * w2.par.T_a * x0[10] / (w2.par.M_G_a * w2.par.V_a);
    parameter Real KcPCA2 = 3e-5;
    parameter Real TiPCA2 = 600;
    parameter Real x0PCA2 = -TiPCA2*u0PCA2/KcPCA2;

    components.PI PCA2(Kc=KcPCA2,Ti=TiPCA2, maxSetPoint=115e5, minSetPoint=75e5, satMax=1, satMin=0, intError(start=x0PCA2, nominal=3.02436e7), u(start=u0PCA2));

    parameter Real u0PCT2(min = 0, max = 1) = (y0PCA2 - PCA2.minSetPoint)  / (PCA2.maxSetPoint - PCA2.minSetPoint);
    parameter Real ZcT2 = w2.par.Z_ct; //x0[9]/100e5; //1 + 1e-5*(x0[9] - w1.par.P_z)*w2.par.Slope_z;
    parameter Real y0PCT2(min = 0) = ZcT2*x0[11]*R*w2.par.T_r/((w2.V_r - (x0[12]-w2.par.rho_L*w2.par.L_bh*w2.par.S_bh)/w2.par.rho_L)*w2.par.M_G_r_t) - y0PCINL - m.Fric;
    Modelica.Blocks.Interfaces.RealInput valveDP_whd2_SetPoint(start=y0PCT2);
    parameter Real KcPCT2 = -4e-9;
    parameter Real TiPCT2 = 100;
    parameter Real x0PCT2 = -TiPCT2*u0PCT2/KcPCT2;

    components.PI PCT2(Kc=KcPCT2,Ti=TiPCT2, satMax=1, satMin=0, intError(start=x0PCT2, nominal=8.86107e9), u(start = u0PCT2));

    //components.firstOrder SPTOPFilter(k=1, T=300, y(start = y0PCTOP, nominal = 1e5));

    //components.firstOrder PCT1Filter(k=1, T=150, y(start = u0PCT1, nominal = 1));
    //components.firstOrder PCT2Filter(k=1, T=150, y(start = u0PCT2, nominal = 1));
    //components.firstOrder PCTOPFilter(k=1, T=150, y(start = u0PCTOP, nominal = 1));

    //components.firstOrder PCA1Filter(k=1, T=0.1, y(start = u0PCA1, nominal = 1));
    //components.firstOrder PCA2Filter(k=1, T=0.1, y(start = u0PCA2, nominal = 1));
    //components.firstOrder PCINLFilter(k=1, T=0.1, y(start = u0PCINL, nominal = 1));

  equation
  connect(valveDP_whd1_SetPoint,PCT1.extSetPoint);
  PCT1.measurement = w1.deltaP_valve;

  connect(valveDP_whd2_SetPoint,PCT2.extSetPoint);
  PCT2.measurement = w2.deltaP_valve;

  //valveDP_top_SetPoint = SPTOPFilter.y;
  connect(valveDP_top_SetPoint,PCTOP.extSetPoint);
  PCTOP.measurement = p.deltaP_valve;

  //connect(PCT1.u,PCT1Filter.u);
  //connect(PCT1Filter.y,PCA1.extSetPoint);
  connect(PCT1.u,PCA1.extSetPoint);
  PCA1.measurement = w1.fin.p;

  //connect(PCT2.u,PCT2Filter.u);
  //connect(PCT2Filter.y,PCA2.extSetPoint);
  connect(PCT2.u,PCA2.extSetPoint);
  PCA2.measurement = w2.fin.p;

  //connect(PCTOP.u,PCTOPFilter.u);
  //connect(PCTOPFilter.y,PCINL.extSetPoint);
  connect(PCTOP.u,PCINL.extSetPoint);
  PCINL.measurement = p.fin.p;

  connect(PCA1.u,w1.z1);
  connect(PCA2.u,w2.z1);
  connect(PCINL.u,p.z);

  //connect(PCA1.u,PCA1Filter.u);
  //connect(PCA2.u,PCA2Filter.u);
  //connect(PCINL.u,PCINLFilter.u);

  //connect(PCA1Filter.y,w1.z1);
  //connect(PCA2Filter.y,w2.z1);
  //connect(PCINL.u,p.z);
  //connect(PCINLFilter.y,p.z);

  //connect(pipelinePressureSetPoint,p.Pin_ss);

    annotation (experiment(
        StopTime=36000,
        Interval=1,
        Tolerance=1e-012,
        Algorithm="Dassl"), __Dymola_experimentSetupOutput);
  end netPres_PIDstruc2;

  model netPres_PIDstruc2_Simulate

    replaceable package components = networkModels.networkComponents;

    import SI = Modelica.SIunits;
    import Modelica.Constants.R;
    import Modelica.Constants.pi;

    networkModels.gasliftNetworks.netPres_PIDstruc2 netPid;

    //parameter Real[3] PresssureSetpoints = { 2.3244e+006, 2.55069e+006, 2.56983e+006};

    //2.54381000e+06;
    parameter SI.Pressure PCT1SetPoint= netPid.ZcT1*netPid.x0[5]*R*netPid.w1.par.T_r/((netPid.w1.V_r - (netPid.x0[6]-netPid.w1.par.rho_L*netPid.w1.par.L_bh*netPid.w1.par.S_bh)/netPid.w1.par.rho_L)*netPid.w1.par.M_G_r_t)
                                        - netPid.y0PCINL - netPid.m.Fric;

    //2.57117000e+06;
    parameter SI.Pressure PCT2SetPoint= netPid.ZcT2*netPid.x0[11]*R*netPid.w2.par.T_r/((netPid.w2.V_r - (netPid.x0[12]-netPid.w2.par.rho_L*netPid.w2.par.L_bh*netPid.w2.par.S_bh)/netPid.w2.par.rho_L)*netPid.w2.par.M_G_r_t)
                                          - netPid.y0PCINL - netPid.m.Fric;

    //2.27915000e+06;
    parameter SI.Pressure PCTOPSetPoint = netPid.ZcTOP*netPid.x0[16]*R*netPid.p.par.T2/(netPid.p.par.M_Gr*((pi * netPid.p.par.r2 ^ 2) * (netPid.p.par.L2 + netPid.p.par.L3) - netPid.x0[17] / netPid.p.par.rho_L)) - netPid.sep.Psep_nom;

     components.gasManifoldStep gm1(offset = 1.28023,  height = 0.07, startTime=2*3600, fout(p(start=86.5e5)));
     components.gasManifoldStep gm2(offset = 1.30082,  height = 0.07, startTime=4*3600, fout(p(start=91.6e5)));
     //Modelica.Blocks.Sources.Step gm1_step(offset = 1,  height = 0.07, startTime=4*3600);
     //Modelica.Blocks.Sources.Step gm2_step(offset = 1,  height = 0.07, startTime=4*3600);

     //Modelica.Blocks.Sources.Constant up(k = PCTOPSetPoint);
     Modelica.Blocks.Sources.Step step1(offset = PCTOPSetPoint,  height = -0.05e5, startTime=6*3600);
     Modelica.Blocks.Sources.Step step2(offset = 0,  height = -0.05e5, startTime=9*3600);

     //components.gasManifold gm_1;
     //components.gasManifold gm_2;

     Modelica.Blocks.Math.Add u3(u1=step1.y, k1=1, u2=step2.y, k2=1);

     Modelica.Blocks.Sources.Step u4(offset = PCT1SetPoint,  height = -0.1e5, startTime=7*3600);
     Modelica.Blocks.Sources.Step u5(offset = PCT2SetPoint,  height = -0.1e5, startTime=8*3600);

     Modelica.Blocks.Sources.Constant GOR1(k = 0);
     Modelica.Blocks.Sources.Constant GOR2(k = 0);

     //Modelica.Blocks.Sources.Constant Pr1(k = 160);
     //Modelica.Blocks.Sources.Constant Pr2(k = 170);

     Modelica.Blocks.Sources.Constant d_PI1(k = 0);
     Modelica.Blocks.Sources.Constant d_PI2(k = 0);

     Modelica.Blocks.Sources.Constant d_Psep(k = 0);

  equation
     connect(GOR1.y, netPid.w1.GOR);
     connect(GOR2.y, netPid.w2.GOR);
     connect(d_PI1.y, netPid.w1.d_PI);
     connect(d_PI2.y, netPid.w2.d_PI);
     //connect(Pr1.y, netPid.w1.P_res);
     //connect(Pr2.y, netPid.w2.P_res);
     connect(gm1.fout, netPid.w1.fin);
     connect(gm2.fout, netPid.w2.fin);
     connect(u3.y, netPid.valveDP_top_SetPoint);
     connect(u4.y, netPid.valveDP_whd1_SetPoint);
     connect(u5.y, netPid.valveDP_whd2_SetPoint);
     connect(d_Psep.y,netPid.sep.d_Psep)
    annotation (experiment(
        StopTime=36000,
        Interval=10,
        Tolerance=1e-008,
        Algorithm="Dassl"), __Dymola_experimentSetupOutput);

    annotation (experiment(
        StopTime=36000,
        Interval=10,
        Tolerance=1e-008,
        Algorithm="Dassl"));
  end netPres_PIDstruc2_Simulate;

  partial model netSim

    import SI = Modelica.SIunits;

    replaceable package components = networkModels.networkComponents;

    parameter Real[15] x0 = {1.29622, 9.9943e+006, 4855.59, 337.873, 8066.91, 1.32441, 1.04289e+007, 5028.53, 330.103, 8549.66, 2.543e+006, 1812.14, 26909, 46.412, 1563.45};

    SI.Pressure Dp_w1(min=0,nominal=1e5);
    SI.Pressure Dp_w2(min=0,nominal=1e5);

    networkModels.networkComponents.gasliftWell_Sim
                               w1sim(par=networkModels.Instances.well1(P_z=x0[2]),
                              w_G_f(start=x0[1], nominal=x0[1]),
                              P_bh_f(start=x0[2], nominal=x0[2]),
                              m_Ga(start=x0[3], nominal=x0[3]),
                              m_Gw(start=x0[4], nominal=x0[4]),
                              m_Lw(start=x0[5], nominal=x0[5]));

    networkModels.networkComponents.gasliftWell_Sim
                               w2sim(par=networkModels.Instances.well2(P_z=x0[7]),
                              w_G_f(start=x0[6], nominal=x0[6]),
                              P_bh_f(start=x0[7], nominal=x0[7]),
                              m_Ga(start=x0[8], nominal=x0[8]),
                              m_Gw(start=x0[9], nominal=x0[9]),
                              m_Lw(start=x0[10], nominal=x0[10]));

    components.manifold m( Fric=18105);

    components.pipelineRiser p(P1_f(start=x0[11], nominal=x0[11]),
                               m_gp(start=x0[12], nominal=x0[12]),
                               m_lp(start=x0[13], nominal=x0[13]),
                               m_gr(start=x0[14], nominal=x0[14]),
                               m_lr(start=x0[15], nominal=x0[15]));

    components.separatorTwoPhase sep(Psep(start=sep.Psep_nom),fin(w(start={2.7,36.3})));

  equation
    Dp_w1 = w1sim.P_r_t - p.fin.p;
    Dp_w2 = w2sim.P_r_t - p.fin.p;

    connect(w1sim.fout, m.fin);
    connect(w2sim.fout, m.fin);
    connect(m.fout, p.fin);
    connect(p.fout, sep.fin);
  end netSim;

  model netSim_PIDstruc2
      extends networkModels.gasliftNetworks.netSim;

    import SI = Modelica.SIunits;
    import Modelica.Constants.R;
    import Modelica.Constants.pi;

    parameter Real u0PCINL(min = 0, max = 1) = 0.577059;
    parameter Real ZcINL = p.par.Z_c0; //x0[13]/p.par.P1; //1 + 1e-5*(x0[13] - p.par.P_z)*p.par.Slope_z;
    parameter Real y0PCINL(min = 2e6, max = 4e6) = ZcINL * x0[12] * R * p.par.T1 / (p.par.M_Gp * (p.V1 - x0[13] / p.par.rho_L));
    //Modelica.Blocks.Interfaces.RealInput pipelinePressureSetPoint(start=y0PCINL);
    parameter Real KcPCINL = 4e-6;
    parameter Real TiPCINL = 120;
    parameter Real x0PCINL = -TiPCINL*u0PCINL/KcPCINL;

    components.PI PCINL(Kc=KcPCINL, Ti=TiPCINL, maxSetPoint=40e5, minSetPoint=20e5, satMax=1,satMin=0,intError(start=x0PCINL, nominal=3.00778e7), u(start = u0PCINL));

    parameter Real u0PCTOP(min = 0, max = 1) = (y0PCINL - PCINL.minSetPoint)  / (PCINL.maxSetPoint - PCINL.minSetPoint);
    parameter Real ZcTOP = p.par.Z_c0; //x0[13]/p.par.P1; //1 + 1e-5*(x0[13] - p.par.P_z)*p.par.Slope_z;
    parameter Real y0PCTOP(min = 0) = ZcTOP*x0[14]*R*p.par.T2/(p.par.M_Gr*((pi * p.par.r2 ^ 2) * (p.par.L2 + p.par.L3) - x0[15] / p.par.rho_L)) - sep.Psep_nom;
    Modelica.Blocks.Interfaces.RealInput valveDP_top_SetPoint(start=y0PCTOP);
    parameter Real KcPCTOP = -2e-7;
    parameter Real TiPCTOP = 300;
    parameter Real x0PCTOP = -TiPCTOP*u0PCTOP/KcPCTOP;

    components.PI PCTOP(Kc=KcPCTOP, Ti=TiPCTOP, satMax=1,satMin=0,intError(start=x0PCTOP, nominal=1.08894e8), u(start = u0PCTOP));

    parameter Real u0PCA1(min = 0, max = 1) = 0.580743;
    parameter Real ZcA1 = w1sim.par.Z_ca; //x0[3]/100e5;//1 + 1e-5*(x0[3] - w1sim.par.P_z)*w1sim.par.Slope_z;
    parameter Real y0PCA1(min = 60e5, max = 100e5) = ZcA1 * R * w1sim.par.T_a * x0[3] / (w1sim.par.M_G_a * w1sim.par.V_a);
    parameter Real KcPCA1 = 3e-5;
    parameter Real TiPCA1 = 600;
    parameter Real x0PCA1 = -TiPCA1*u0PCA1/KcPCA1;

    components.PI PCA1(Kc=KcPCA1,Ti=TiPCA1, maxSetPoint=110e5, minSetPoint=70e5, satMax=1, satMin=0, intError(start=x0PCA1, nominal=3.02824e7), u(start=u0PCA1));

    parameter Real u0PCT1(min = 0, max = 1) = (y0PCA1 - PCA1.minSetPoint)  / (PCA1.maxSetPoint - PCA1.minSetPoint);
    parameter Real ZcT1 = w1sim.par.Z_ct; //x0[3]/100e5; //1 + 1e-5*(x0[3] - w1sim.par.P_z)*w1sim.par.Slope_z;
    parameter Real y0PCT1(min = 0) = ZcT1*x0[4]*R*w1sim.par.T_r/((w1sim.V_r - (x0[5]-w1sim.par.rho_L*w1sim.par.L_bh*w1sim.par.S_bh)/w1sim.par.rho_L)*w1sim.par.M_G_r_t) - y0PCINL - m.Fric;
    Modelica.Blocks.Interfaces.RealInput valveDP_whd1_SetPoint(start=y0PCT1);
    parameter Real KcPCT1 = -4e-9;
    parameter Real TiPCT1 = 100;
    parameter Real x0PCT1 = -TiPCT1*u0PCT1/KcPCT1;

    components.PI PCT1(Kc=KcPCT1,Ti=TiPCT1, satMax=1,satMin=0,intError(start=x0PCT1, nominal=9.25742e9), u(start = u0PCT1));

    parameter Real u0PCA2(min = 0, max = 1) = 0.608612;
    parameter Real ZcA2 = w2sim.par.Z_ca; //x0[9]/100e5; //1 + 1e-5*(x0[9] - w1sim.par.P_z)*w2sim.par.Slope_z;
    parameter Real y0PCA2(min = 60e5, max = 100e5) = ZcA2 * R * w2sim.par.T_a * x0[8] / (w2sim.par.M_G_a * w2sim.par.V_a);
    parameter Real KcPCA2 = 3e-5;
    parameter Real TiPCA2 = 600;
    parameter Real x0PCA2 = -TiPCA2*u0PCA2/KcPCA2;

    components.PI PCA2(Kc=KcPCA2,Ti=TiPCA2, maxSetPoint=115e5, minSetPoint=75e5, satMax=1, satMin=0, intError(start=x0PCA2, nominal=3.02436e7), u(start=u0PCA2));

    parameter Real u0PCT2(min = 0, max = 1) = (y0PCA2 - PCA2.minSetPoint)  / (PCA2.maxSetPoint - PCA2.minSetPoint);
    parameter Real ZcT2 = w2sim.par.Z_ct; //x0[9]/100e5; //1 + 1e-5*(x0[9] - w1sim.par.P_z)*w2sim.par.Slope_z;
    parameter Real y0PCT2(min = 0) = ZcT2*x0[9]*R*w2sim.par.T_r/((w2sim.V_r - (x0[10]-w2sim.par.rho_L*w2sim.par.L_bh*w2sim.par.S_bh)/w2sim.par.rho_L)*w2sim.par.M_G_r_t) - y0PCINL - m.Fric;
    Modelica.Blocks.Interfaces.RealInput valveDP_whd2_SetPoint(start=y0PCT2);
    parameter Real KcPCT2 = -4e-9;
    parameter Real TiPCT2 = 100;
    parameter Real x0PCT2 = -TiPCT2*u0PCT2/KcPCT2;

    components.PI PCT2(Kc=KcPCT2,Ti=TiPCT2, satMax=1, satMin=0, intError(start=x0PCT2, nominal=8.86107e9), u(start = u0PCT2));

  equation
    connect(valveDP_whd1_SetPoint,PCT1.extSetPoint);
    PCT1.measurement = w1sim.deltaP_valve;

    connect(valveDP_whd2_SetPoint,PCT2.extSetPoint);
    PCT2.measurement = w2sim.deltaP_valve;

    connect(valveDP_top_SetPoint,PCTOP.extSetPoint);
    PCTOP.measurement = p.deltaP_valve;

    connect(PCT1.u,PCA1.extSetPoint);
    PCA1.measurement = w1sim.fin.p;

    connect(PCT2.u,PCA2.extSetPoint);
    PCA2.measurement = w2sim.fin.p;

    connect(PCTOP.u,PCINL.extSetPoint);
    PCINL.measurement = p.fin.p;

    connect(PCA1.u,w1sim.z1);
    connect(PCA2.u,w2sim.z1);
    connect(PCINL.u,p.z);

    //connect(pipelinePressureSetPoint,p.Pin_ss);

    annotation (experiment(
        StopTime=36000,
        Interval=1,
        Tolerance=1e-012,
        Algorithm="Dassl"), __Dymola_experimentSetupOutput);
  end netSim_PIDstruc2;

  model netSim_PIDstruc2_Simulate

    replaceable package components = networkModels.networkComponents;

    import SI = Modelica.SIunits;
    import Modelica.Constants.R;
    import Modelica.Constants.pi;

    networkModels.gasliftNetworks.netSim_PIDstruc2 netPid;

    //parameter Real[3] PresssureSetpoints = { 2.3244e+006, 2.55069e+006, 2.56983e+006};

    //2.54381000e+06;
    parameter SI.Pressure PCT1SetPoint= netPid.ZcT1*netPid.x0[4]*R*netPid.w1sim.par.T_r/((netPid.w1sim.V_r - (netPid.x0[5]-netPid.w1sim.par.rho_L*netPid.w1sim.par.L_bh*netPid.w1sim.par.S_bh)/netPid.w1sim.par.rho_L)*netPid.w1sim.par.M_G_r_t)
                                        - netPid.y0PCINL - netPid.m.Fric;

    //2.57117000e+06;
    parameter SI.Pressure PCT2SetPoint= netPid.ZcT2*netPid.x0[9]*R*netPid.w2sim.par.T_r/((netPid.w2sim.V_r - (netPid.x0[10]-netPid.w2sim.par.rho_L*netPid.w2sim.par.L_bh*netPid.w2sim.par.S_bh)/netPid.w2sim.par.rho_L)*netPid.w2sim.par.M_G_r_t)
                                          - netPid.y0PCINL - netPid.m.Fric;

    //2.27915000e+06;
    parameter SI.Pressure PCTOPSetPoint = netPid.ZcTOP*netPid.x0[14]*R*netPid.p.par.T2/(netPid.p.par.M_Gr*((pi * netPid.p.par.r2 ^ 2) * (netPid.p.par.L2 + netPid.p.par.L3) - netPid.x0[15] / netPid.p.par.rho_L)) - netPid.sep.Psep_nom;
     components.gasManifoldStep gm1(offset = 1.29633,  height = 0.03, startTime=2*3600, fout(p(start=86.5e5)));
     components.gasManifoldStep gm2(offset = 1.32452,  height = 0.03, startTime=4*3600, fout(p(start=91.6e5)));
     //Modelica.Blocks.Sources.Step gm1_step(offset = 1,  height = 0.07, startTime=4*3600);
     //Modelica.Blocks.Sources.Step gm2_step(offset = 1,  height = 0.07, startTime=4*3600);

     //Modelica.Blocks.Sources.Constant up(k = PCTOPSetPoint);
     Modelica.Blocks.Sources.Step step1(offset = PCTOPSetPoint,  height = -0.05e5, startTime=6*3600);
     Modelica.Blocks.Sources.Step step2(offset = 0,  height = -0.05e5, startTime=9*3600);

     //components.gasManifold gm_1;
     //components.gasManifold gm_2;

     Modelica.Blocks.Math.Add u3(u1=step1.y, k1=1, u2=step2.y, k2=1);

     Modelica.Blocks.Sources.Step u4(offset = PCT1SetPoint,  height = -0.1e5, startTime=7*3600);
     Modelica.Blocks.Sources.Step u5(offset = PCT2SetPoint,  height = -0.1e5, startTime=8*3600);

     //Modelica.Blocks.Sources.Constant GOR1(k = 0);
     //Modelica.Blocks.Sources.Constant GOR2(k = 0);

     Modelica.Blocks.Sources.Step GOR1(offset = 0, height = 0.005, startTime=2*3600);
     Modelica.Blocks.Sources.Step GOR2(offset = 0, height = 0.005, startTime=4*3600);

     //Modelica.Blocks.Sources.Constant Pr1(k = 160);
     //Modelica.Blocks.Sources.Constant Pr2(k = 170);
     Modelica.Blocks.Sources.Step Pr1(offset = 160,  height = -1, startTime=5*3600);
     Modelica.Blocks.Sources.Step Pr2(offset = 170,  height = -1, startTime=9*3600);

     Modelica.Blocks.Sources.Constant d_PI1(k = 0);
     Modelica.Blocks.Sources.Constant d_PI2(k = 0);

     Modelica.Blocks.Sources.Constant d_Psep(k = 0);

  equation
     connect(GOR1.y, netPid.w1sim.alpha_Gm);
     connect(GOR2.y, netPid.w2sim.alpha_Gm);
     connect(d_PI1.y, netPid.w1sim.d_PI);
     connect(d_PI2.y, netPid.w2sim.d_PI);
     connect(Pr1.y, netPid.w1sim.P_res);
     connect(Pr2.y, netPid.w2sim.P_res);
     connect(gm1.fout, netPid.w1sim.fin);
     connect(gm2.fout, netPid.w2sim.fin);
     connect(u3.y, netPid.valveDP_top_SetPoint);
     connect(u4.y, netPid.valveDP_whd1_SetPoint);
     connect(u5.y, netPid.valveDP_whd2_SetPoint);
     connect(d_Psep.y,netPid.sep.d_Psep)
  annotation (experiment(
      StopTime=36000,
      Interval=10,
      Tolerance=1e-008,
      Algorithm="Dassl"));

  end netSim_PIDstruc2_Simulate;

  partial model netGOR

    import SI = Modelica.SIunits;

    replaceable package components = networkModels.networkComponents;

    //parameter Real[17] x0 = {0, 1.28004, 9.99901e+006, 4745.12, 332.637, 8060.01, 0, 1.30071, 1.04357e+007, 4917.94, 323.47, 8538.99, 2.51954e+006, 1795.42, 26909.2, 46.339, 1576.8};

    parameter Real[17] x0 = {0, 1.29633, 9.99415e+006, 4855.59, 337.87, 8066.59, 0, 1.32452, 1.04287e+007, 5028.51, 330.09, 8549.3, 2.54313e+006, 1812.23, 26909, 46.4101, 1563.37};

    SI.Pressure Dp_w1(min=0,nominal=1e5);
    SI.Pressure Dp_w2(min=0,nominal=1e5);

    networkModels.networkComponents.gasliftWell_GOR
                           w1(par=networkModels.Instances.well1(P_z=x0[3]),
                              alpha_G_m_in(start=x0[1], nominal=0.04),
                              w_G_f(start=x0[2], nominal=x0[2]),
                              P_bh_f(start=x0[3], nominal=x0[3]),
                              m_Ga(start=x0[4], nominal=x0[4]),
                              m_Gw(start=x0[5], nominal=x0[5]),
                              m_Lw(start=x0[6], nominal=x0[6]));

    networkModels.networkComponents.gasliftWell_GOR
                           w2(par=networkModels.Instances.well2(P_z=x0[9]),
                              alpha_G_m_in(start=x0[7], nominal=0.04),
                              w_G_f(start=x0[8], nominal=x0[8]),
                              P_bh_f(start=x0[9], nominal=x0[9]),
                              m_Ga(start=x0[10], nominal=x0[10]),
                              m_Gw(start=x0[11], nominal=x0[11]),
                              m_Lw(start=x0[12], nominal=x0[12]));

    components.manifold m( Fric=18105);

    components.pipelineRiser p(P1_f(start=x0[13], nominal=x0[13]),
                               m_gp(start=x0[14], nominal=x0[14]),
                               m_lp(start=x0[15], nominal=x0[15]),
                               m_gr(start=x0[16], nominal=x0[16]),
                               m_lr(start=x0[17], nominal=x0[17]));

    components.separatorTwoPhase sep(Psep(start=sep.Psep_nom),fin(w(start={2.7,36.3})));

  equation
    Dp_w1 = w1.P_r_t - p.fin.p;
    Dp_w2 = w2.P_r_t - p.fin.p;

    connect(w1.fout, m.fin);
    connect(w2.fout, m.fin);
    connect(m.fout, p.fin);
    connect(p.fout, sep.fin);
  end netGOR;

  model netGOR_PIDstruc2
    extends networkModels.gasliftNetworks.netGOR;

    import SI = Modelica.SIunits;
    import Modelica.Constants.R;
    import Modelica.Constants.pi;

    parameter Real u0PCINL(min = 0, max = 1) = 0.577059;
    parameter Real ZcINL = p.par.Z_c0; //x0[13]/p.par.P1; //1 + 1e-5*(x0[13] - p.par.P_z)*p.par.Slope_z;
    parameter Real y0PCINL(min = 2e6, max = 4e6) = ZcINL * x0[14] * R * p.par.T1 / (p.par.M_Gp * (p.V1 - x0[15] / p.par.rho_L));
    //Modelica.Blocks.Interfaces.RealInput pipelinePressureSetPoint(start=y0PCINL);
    parameter Real KcPCINL = 4e-6;
    parameter Real TiPCINL = 120;
    parameter Real x0PCINL = -TiPCINL*u0PCINL/KcPCINL;

    components.PI PCINL(Kc=KcPCINL, Ti=TiPCINL, maxSetPoint=40e5, minSetPoint=20e5, satMax=1,satMin=0,intError(start=x0PCINL, nominal=3.00778e7), u(start = u0PCINL));

    parameter Real u0PCTOP(min = 0, max = 1) = (y0PCINL - PCINL.minSetPoint)  / (PCINL.maxSetPoint - PCINL.minSetPoint);
    parameter Real ZcTOP = p.par.Z_c0; //x0[13]/p.par.P1; //1 + 1e-5*(x0[13] - p.par.P_z)*p.par.Slope_z;
    parameter Real y0PCTOP(min = 0) = ZcTOP*x0[16]*R*p.par.T2/(p.par.M_Gr*((pi * p.par.r2 ^ 2) * (p.par.L2 + p.par.L3) - x0[17] / p.par.rho_L)) - sep.Psep_nom;
    Modelica.Blocks.Interfaces.RealInput valveDP_top_SetPoint(start=y0PCTOP);
    parameter Real KcPCTOP = -2e-7;
    parameter Real TiPCTOP = 300;
    parameter Real x0PCTOP = -TiPCTOP*u0PCTOP/KcPCTOP;

    components.PI PCTOP(Kc=KcPCTOP, Ti=TiPCTOP, satMax=1,satMin=0,intError(start=x0PCTOP, nominal=1.08894e8), u(start = u0PCTOP));

    parameter Real u0PCA1(min = 0, max = 1) = 0.580743;
    parameter Real ZcA1 = w1.par.Z_ca; //x0[3]/100e5;//1 + 1e-5*(x0[3] - w1.par.P_z)*w1.par.Slope_z;
    parameter Real y0PCA1(min = 60e5, max = 100e5) = ZcA1 * R * w1.par.T_a * x0[4] / (w1.par.M_G_a * w1.par.V_a);
    parameter Real KcPCA1 = 3e-5;
    parameter Real TiPCA1 = 600;
    parameter Real x0PCA1 = -TiPCA1*u0PCA1/KcPCA1;

    components.PI PCA1(Kc=KcPCA1,Ti=TiPCA1, maxSetPoint=110e5, minSetPoint=70e5, satMax=1, satMin=0, intError(start=x0PCA1, nominal=3.02824e7), u(start=u0PCA1));

    parameter Real u0PCT1(min = 0, max = 1) = (y0PCA1 - PCA1.minSetPoint)  / (PCA1.maxSetPoint - PCA1.minSetPoint);
    parameter Real ZcT1 = w1.par.Z_ct; //x0[3]/100e5; //1 + 1e-5*(x0[3] - w1.par.P_z)*w1.par.Slope_z;
    parameter Real y0PCT1(min = 0) = ZcT1*x0[5]*R*w1.par.T_r/((w1.V_r - (x0[6]-w1.par.rho_L*w1.par.L_bh*w1.par.S_bh)/w1.par.rho_L)*w1.par.M_G_r_t) - y0PCINL - m.Fric;
    Modelica.Blocks.Interfaces.RealInput valveDP_whd1_SetPoint(start=y0PCT1);
    parameter Real KcPCT1 = -4e-9;
    parameter Real TiPCT1 = 100;
    parameter Real x0PCT1 = -TiPCT1*u0PCT1/KcPCT1;

    components.PI PCT1(Kc=KcPCT1,Ti=TiPCT1, satMax=1,satMin=0,intError(start=x0PCT1, nominal=9.25742e9), u(start = u0PCT1));

    parameter Real u0PCA2(min = 0, max = 1) = 0.608612;
    parameter Real ZcA2 = w2.par.Z_ca; //x0[9]/100e5; //1 + 1e-5*(x0[9] - w1.par.P_z)*w2.par.Slope_z;
    parameter Real y0PCA2(min = 60e5, max = 100e5) = ZcA2 * R * w2.par.T_a * x0[10] / (w2.par.M_G_a * w2.par.V_a);
    parameter Real KcPCA2 = 3e-5;
    parameter Real TiPCA2 = 600;
    parameter Real x0PCA2 = -TiPCA2*u0PCA2/KcPCA2;

    components.PI PCA2(Kc=KcPCA2,Ti=TiPCA2, maxSetPoint=115e5, minSetPoint=75e5, satMax=1, satMin=0, intError(start=x0PCA2, nominal=3.02436e7), u(start=u0PCA2));

    parameter Real u0PCT2(min = 0, max = 1) = (y0PCA2 - PCA2.minSetPoint)  / (PCA2.maxSetPoint - PCA2.minSetPoint);
    parameter Real ZcT2 = w2.par.Z_ct; //x0[9]/100e5; //1 + 1e-5*(x0[9] - w1.par.P_z)*w2.par.Slope_z;
    parameter Real y0PCT2(min = 0) = ZcT2*x0[11]*R*w2.par.T_r/((w2.V_r - (x0[12]-w2.par.rho_L*w2.par.L_bh*w2.par.S_bh)/w2.par.rho_L)*w2.par.M_G_r_t) - y0PCINL - m.Fric;
    Modelica.Blocks.Interfaces.RealInput valveDP_whd2_SetPoint(start=y0PCT2);
    parameter Real KcPCT2 = -4e-9;
    parameter Real TiPCT2 = 100;
    parameter Real x0PCT2 = -TiPCT2*u0PCT2/KcPCT2;

    components.PI PCT2(Kc=KcPCT2,Ti=TiPCT2, satMax=1, satMin=0, intError(start=x0PCT2, nominal=8.86107e9), u(start = u0PCT2));

    //components.firstOrder SPTOPFilter(k=1, T=300, y(start = y0PCTOP, nominal = 1e5));

    //components.firstOrder PCT1Filter(k=1, T=30, y(start = u0PCT1, nominal = 1));
    //components.firstOrder PCT2Filter(k=1, T=30, y(start = u0PCT2, nominal = 1));
    //components.firstOrder PCTOPFilter(k=1, T=30, y(start = u0PCTOP, nominal = 1));

    //components.firstOrder PCA1Filter(k=1, T=0.1, y(start = u0PCA1, nominal = 1));
    //components.firstOrder PCA2Filter(k=1, T=0.1, y(start = u0PCA2, nominal = 1));
    //components.firstOrder PCINLFilter(k=1, T=0.1, y(start = u0PCINL, nominal = 1));

  equation
  connect(valveDP_whd1_SetPoint,PCT1.extSetPoint);
  PCT1.measurement = w1.deltaP_valve;

  connect(valveDP_whd2_SetPoint,PCT2.extSetPoint);
  PCT2.measurement = w2.deltaP_valve;

  //valveDP_top_SetPoint = SPTOPFilter.y;
  connect(valveDP_top_SetPoint,PCTOP.extSetPoint);
  PCTOP.measurement = p.deltaP_valve;

  //connect(PCT1.u,PCT1Filter.u);
  //connect(PCT1Filter.y,PCA1.extSetPoint);
  connect(PCT1.u,PCA1.extSetPoint);
  PCA1.measurement = w1.fin.p;

  //connect(PCT2.u,PCT2Filter.u);
  //connect(PCT2Filter.y,PCA2.extSetPoint);
  connect(PCT2.u,PCA2.extSetPoint);
  PCA2.measurement = w2.fin.p;

  //connect(PCTOP.u,PCTOPFilter.u);
  //connect(PCTOPFilter.y,PCINL.extSetPoint);
  connect(PCTOP.u,PCINL.extSetPoint);
  PCINL.measurement = p.fin.p;

  connect(PCA1.u,w1.z1);
  connect(PCA2.u,w2.z1);
  connect(PCINL.u,p.z);

  //connect(PCA1.u,PCA1Filter.u);
  //connect(PCA2.u,PCA2Filter.u);
  //connect(PCINL.u,PCINLFilter.u);

  //connect(PCA1Filter.y,w1.z1);
  //connect(PCA2Filter.y,w2.z1);
  //connect(PCINL.u,p.z);
  //connect(PCINLFilter.y,p.z);

  //connect(pipelinePressureSetPoint,p.Pin_ss);

    annotation (experiment(
        StopTime=36000,
        Interval=1,
        Tolerance=1e-012,
        Algorithm="Dassl"), __Dymola_experimentSetupOutput);
  end netGOR_PIDstruc2;

  model netGOR_PIDstruc2_Simulate

    replaceable package components = networkModels.networkComponents;

    import SI = Modelica.SIunits;
    import Modelica.Constants.R;
    import Modelica.Constants.pi;

    networkModels.gasliftNetworks.netGOR_PIDstruc2 netPid;

    //parameter Real[3] PresssureSetpoints = { 2.3244e+006, 2.55069e+006, 2.56983e+006};

    //2.54381000e+06;
  parameter SI.Pressure PCT1SetPoint= netPid.ZcT1*netPid.x0[5]*R*netPid.w1.par.T_r/((netPid.w1.V_r - (netPid.x0[6]-netPid.w1.par.rho_L*netPid.w1.par.L_bh*netPid.w1.par.S_bh)/netPid.w1.par.rho_L)*netPid.w1.par.M_G_r_t)
                                        - netPid.y0PCINL - netPid.m.Fric;

    //2.57117000e+06;
    parameter SI.Pressure PCT2SetPoint= netPid.ZcT2*netPid.x0[11]*R*netPid.w2.par.T_r/((netPid.w2.V_r - (netPid.x0[12]-netPid.w2.par.rho_L*netPid.w2.par.L_bh*netPid.w2.par.S_bh)/netPid.w2.par.rho_L)*netPid.w2.par.M_G_r_t)
                                          - netPid.y0PCINL - netPid.m.Fric;

    //2.27915000e+06;
    parameter SI.Pressure PCTOPSetPoint = netPid.ZcTOP*netPid.x0[16]*R*netPid.p.par.T2/(netPid.p.par.M_Gr*((pi * netPid.p.par.r2 ^ 2) * (netPid.p.par.L2 + netPid.p.par.L3) - netPid.x0[17] / netPid.p.par.rho_L)) - netPid.sep.Psep_nom;

     components.gasManifoldStep gm1(offset = 1.29633,  height = 0.03, startTime=2*3600, fout(p(start=86.5e5)));
     components.gasManifoldStep gm2(offset = 1.32452,  height = 0.03, startTime=4*3600, fout(p(start=91.6e5)));
     //Modelica.Blocks.Sources.Step gm1_step(offset = 1,  height = 0.07, startTime=4*3600);
     //Modelica.Blocks.Sources.Step gm2_step(offset = 1,  height = 0.07, startTime=4*3600);

     //Modelica.Blocks.Sources.Constant up(k = PCTOPSetPoint);
     Modelica.Blocks.Sources.Step step1(offset = PCTOPSetPoint,  height = -0.05e5, startTime=6*3600);
     Modelica.Blocks.Sources.Step step2(offset = 0,  height = -0.05e5, startTime=9*3600);

     //components.gasManifold gm_1;
     //components.gasManifold gm_2;

     Modelica.Blocks.Math.Add u3(u1=step1.y, k1=1, u2=step2.y, k2=1);

     Modelica.Blocks.Sources.Step u4(offset = PCT1SetPoint,  height = -0.1e5, startTime=7*3600);
     Modelica.Blocks.Sources.Step u5(offset = PCT2SetPoint,  height = -0.1e5, startTime=8*3600);

     //Modelica.Blocks.Sources.Constant GOR1(k = 0);
     //Modelica.Blocks.Sources.Constant GOR2(k = 0);

     Modelica.Blocks.Sources.Step GOR1(offset = 0, height = 0.1, startTime=2*3600);
     Modelica.Blocks.Sources.Step GOR2(offset = 0, height = 0.1, startTime=4*3600);

     //Modelica.Blocks.Sources.Constant Pr1(k = 160);
     //Modelica.Blocks.Sources.Constant Pr2(k = 170);

     Modelica.Blocks.Sources.Constant d_PI1(k = 0);
     Modelica.Blocks.Sources.Constant d_PI2(k = 0);

     Modelica.Blocks.Sources.Constant d_Psep(k = 0);

  equation
     connect(GOR1.y, netPid.w1.GOR);
     connect(GOR2.y, netPid.w2.GOR);
     connect(d_PI1.y, netPid.w1.d_PI);
     connect(d_PI2.y, netPid.w2.d_PI);
     //connect(Pr1.y, netPid.w1.P_res);
     //connect(Pr2.y, netPid.w2.P_res);
     connect(gm1.fout, netPid.w1.fin);
     connect(gm2.fout, netPid.w2.fin);
     connect(u3.y, netPid.valveDP_top_SetPoint);
     connect(u4.y, netPid.valveDP_whd1_SetPoint);
     connect(u5.y, netPid.valveDP_whd2_SetPoint);
     connect(d_Psep.y,netPid.sep.d_Psep)
    annotation (experiment(
        StopTime=36000,
        Interval=10,
        Tolerance=1e-008,
        Algorithm="Dassl"), __Dymola_experimentSetupOutput);

    annotation (experiment(
        StopTime=36000,
        Interval=10,
        Tolerance=1e-008,
        Algorithm="Dassl"));
  end netGOR_PIDstruc2_Simulate;

  partial model net

    import SI = Modelica.SIunits;

    replaceable package components = networkModels.networkComponents;

    parameter Real[15] x0 = {1, 1.018e+007, 4763.5, 281.602, 8895.96, 1, 1.0664e+007, 4947.4, 275.389, 9489.86, 2.18077e+006, 1553.96, 26912, 45.5116, 1799.56};

    // {160, 1, 1.0165e+007, 5078.23, 273.556, 7973.3, 170, 1, 1.06757e+007, 5230.69, 267.252, 8657.19, 2.19496e+006, 1411.81, 26432.9, 38.3596, 2130.84};
    // x0 = {160, 1, 1.01845e+007, 5107.89, 274.462, 8208.55, 170, 1, 1.0658e+007, 5241.01, 265.91, 8831.32, 2.08477e+006, 1485.67, 26310.4, 43.4316, 1716.27};

    SI.Pressure Dp_w1(min=0,nominal=1e5);
    SI.Pressure Dp_w2(min=0,nominal=1e5);

    networkModels.networkComponents.gasliftWell
                               w1(par=networkModels.Instances.well1(P_z=9.99901e+006),
                              w_G_f(start=x0[1], nominal=x0[1]),
                              P_bh_f(start=x0[2], nominal=x0[2]),
                              m_Ga(start=x0[3], nominal=x0[3]),
                              m_Gw(start=x0[4], nominal=x0[4]),
                              m_Lw(start=x0[5], nominal=x0[5]));

    networkModels.networkComponents.gasliftWell
                               w2(par=networkModels.Instances.well2(P_z=1.04357e+007),
                              w_G_f(start=x0[6], nominal=x0[6]),
                              P_bh_f(start=x0[7], nominal=x0[7]),
                              m_Ga(start=x0[8], nominal=x0[8]),
                              m_Gw(start=x0[9], nominal=x0[9]),
                              m_Lw(start=x0[10], nominal=x0[10]));

    components.manifold m( Fric=18105);

    components.pipelineRiser p(P1_f(start=x0[11], nominal=x0[11]),
                               m_gp(start=x0[12], nominal=x0[12]),
                               m_lp(start=x0[13], nominal=x0[13]),
                               m_gr(start=x0[14], nominal=x0[14]),
                               m_lr(start=x0[15], nominal=x0[15]));

    components.separatorTwoPhase sep(Psep(start=sep.Psep_nom),fin(w(start={2.7,36.3})));

  equation
    Dp_w1 = w1.P_r_t - p.fin.p;
    Dp_w2 = w2.P_r_t - p.fin.p;

    connect(w1.fout, m.fin);
    connect(w2.fout, m.fin);
    connect(m.fout, p.fin);
    connect(p.fout, sep.fin);
  end net;

  model net_PIDstruc2
    extends networkModels.gasliftNetworks.net;

    import SI = Modelica.SIunits;
    import Modelica.Constants.R;
    import Modelica.Constants.pi;

    parameter Real u0PCINL(min = 0, max = 1) = 0.5;
    parameter Real ZcINL = p.par.Z_c0; //x0[13]/p.par.P1; //1 + 1e-5*(x0[13] - p.par.P_z)*p.par.Slope_z;
    parameter Real y0PCINL(min = 2e6, max = 4e6) = ZcINL * x0[12] * R * p.par.T1 / (p.par.M_Gp * (p.V1 - x0[13] / p.par.rho_L));
    //Modelica.Blocks.Interfaces.RealInput pipelinePressureSetPoint(start=y0PCINL);
    parameter Real KcPCINL = 2e-6;
    parameter Real TiPCINL = 300;
    parameter Real x0PCINL = -TiPCINL*u0PCINL/KcPCINL;

    components.PI PCINL(Kc=KcPCINL, Ti=TiPCINL, maxSetPoint=40e5, minSetPoint=20e5, satMax=1,satMin=0,intError(start=x0PCINL, nominal=3.00778e7), u(start = u0PCINL));

    parameter Real u0PCTOP(min = 0, max = 1) = (y0PCINL - PCINL.minSetPoint)  / (PCINL.maxSetPoint - PCINL.minSetPoint);
    parameter Real ZcTOP = p.par.Z_c0; //x0[13]/p.par.P1; //1 + 1e-5*(x0[13] - p.par.P_z)*p.par.Slope_z;
    parameter Real y0PCTOP(min = 0) = ZcTOP*x0[14]*R*p.par.T2/(p.par.M_Gr*((pi * p.par.r2 ^ 2) * (p.par.L2 + p.par.L3) - x0[15] / p.par.rho_L)) - sep.Psep_nom;
    Modelica.Blocks.Interfaces.RealInput valveDP_top_SetPoint(start=y0PCTOP);
    parameter Real KcPCTOP = -2e-7;
    parameter Real TiPCTOP = 300;
    parameter Real x0PCTOP = -TiPCTOP*u0PCTOP/KcPCTOP;

    components.PI PCTOP(Kc=KcPCTOP, Ti=TiPCTOP, satMax=1,satMin=0,intError(start=x0PCTOP, nominal=1.08894e8), u(start = u0PCTOP));

    parameter Real u0PCA1(min = 0, max = 1) = 0.5;
    parameter Real ZcA1 = w1.par.Z_ca; //x0[3]/100e5;//1 + 1e-5*(x0[3] - w1.par.P_z)*w1.par.Slope_z;
    parameter Real y0PCA1(min = 60e5, max = 100e5) = ZcA1 * R * w1.par.T_a * x0[3] / (w1.par.M_G_a * w1.par.V_a);
    parameter Real KcPCA1 = 2e-5;
    parameter Real TiPCA1 = 600;
    parameter Real x0PCA1 = -TiPCA1*u0PCA1/KcPCA1;

    components.PI PCA1(Kc=KcPCA1,Ti=TiPCA1, maxSetPoint=110e5, minSetPoint=70e5, satMax=1, satMin=0, intError(start=x0PCA1, nominal=3.02824e7), u(start=u0PCA1));

    parameter Real u0PCT1(min = 0, max = 1) = (y0PCA1 - PCA1.minSetPoint)  / (PCA1.maxSetPoint - PCA1.minSetPoint);
    parameter Real ZcT1 = w1.par.Z_ct; //x0[3]/100e5; //1 + 1e-5*(x0[3] - w1.par.P_z)*w1.par.Slope_z;
    parameter Real y0PCT1(min = 0) = ZcT1*x0[4]*R*w1.par.T_r/((w1.V_r - (x0[5]-w1.par.rho_L*w1.par.L_bh*w1.par.S_bh)/w1.par.rho_L)*w1.par.M_G_r_t) - y0PCINL - m.Fric;
    Modelica.Blocks.Interfaces.RealInput valveDP_whd1_SetPoint(start=y0PCT1);
    parameter Real KcPCT1 = -4e-9;
    parameter Real TiPCT1 = 100;
    parameter Real x0PCT1 = -TiPCT1*u0PCT1/KcPCT1;

    components.PI PCT1(Kc=KcPCT1,Ti=TiPCT1, satMax=1,satMin=0,intError(start=x0PCT1, nominal=9.25742e9), u(start = u0PCT1));

    parameter Real u0PCA2(min = 0, max = 1) = 0.5;
    parameter Real ZcA2 = w2.par.Z_ca; //x0[9]/100e5; //1 + 1e-5*(x0[9] - w1.par.P_z)*w2.par.Slope_z;
    parameter Real y0PCA2(min = 60e5, max = 100e5) = ZcA2 * R * w2.par.T_a * x0[8] / (w2.par.M_G_a * w2.par.V_a);
    parameter Real KcPCA2 = 2e-5;
    parameter Real TiPCA2 = 600;
    parameter Real x0PCA2 = -TiPCA2*u0PCA2/KcPCA2;

    components.PI PCA2(Kc=KcPCA2,Ti=TiPCA2, maxSetPoint=115e5, minSetPoint=75e5, satMax=1, satMin=0, intError(start=x0PCA2, nominal=3.02436e7), u(start=u0PCA2));

    parameter Real u0PCT2(min = 0, max = 1) = (y0PCA2 - PCA2.minSetPoint)  / (PCA2.maxSetPoint - PCA2.minSetPoint);
    parameter Real ZcT2 = w2.par.Z_ct; //x0[9]/100e5; //1 + 1e-5*(x0[9] - w1.par.P_z)*w2.par.Slope_z;
    parameter Real y0PCT2(min = 0) = ZcT2*x0[9]*R*w2.par.T_r/((w2.V_r - (x0[10]-w2.par.rho_L*w2.par.L_bh*w2.par.S_bh)/w2.par.rho_L)*w2.par.M_G_r_t) - y0PCINL - m.Fric;
    Modelica.Blocks.Interfaces.RealInput valveDP_whd2_SetPoint(start=y0PCT2);
    parameter Real KcPCT2 = -4e-9;
    parameter Real TiPCT2 = 100;
    parameter Real x0PCT2 = -TiPCT2*u0PCT2/KcPCT2;

    components.PI PCT2(Kc=KcPCT2,Ti=TiPCT2, satMax=1, satMin=0, intError(start=x0PCT2, nominal=8.86107e9), u(start = u0PCT2));

    //components.firstOrder SPTOPFilter(k=1, T=300, y(start = y0PCTOP, nominal = 1e5));

    //components.firstOrder PCT1Filter(k=1, T=150, y(start = u0PCT1, nominal = 1));
    //components.firstOrder PCT2Filter(k=1, T=150, y(start = u0PCT2, nominal = 1));
    //components.firstOrder PCTOPFilter(k=1, T=30, y(start = u0PCTOP, nominal = 1));

    //components.firstOrder PCA1Filter(k=1, T=0.1, y(start = u0PCA1, nominal = 1));
    //components.firstOrder PCA2Filter(k=1, T=0.1, y(start = u0PCA2, nominal = 1));
    //components.firstOrder PCINLFilter(k=1, T=0.1, y(start = u0PCINL, nominal = 1));

  equation
  connect(valveDP_whd1_SetPoint,PCT1.extSetPoint);
  PCT1.measurement = w1.deltaP_valve;

  connect(valveDP_whd2_SetPoint,PCT2.extSetPoint);
  PCT2.measurement = w2.deltaP_valve;

  //valveDP_top_SetPoint = SPTOPFilter.y;
  connect(valveDP_top_SetPoint,PCTOP.extSetPoint);
  PCTOP.measurement = p.deltaP_valve;

  //connect(PCT1.u,PCT1Filter.u);
  //connect(PCT1Filter.y,PCA1.extSetPoint);
  connect(PCT1.u,PCA1.extSetPoint);
  PCA1.measurement = w1.fin.p;

  //connect(PCT2.u,PCT2Filter.u);
  //connect(PCT2Filter.y,PCA2.extSetPoint);
  connect(PCT2.u,PCA2.extSetPoint);
  PCA2.measurement = w2.fin.p;

  //connect(PCTOP.u,PCTOPFilter.u);
  //connect(PCTOPFilter.y,PCINL.extSetPoint);
  connect(PCTOP.u,PCINL.extSetPoint);
  PCINL.measurement = p.fin.p;

  connect(PCA1.u,w1.z1);
  connect(PCA2.u,w2.z1);
  connect(PCINL.u,p.z);

  //connect(PCA1.u,PCA1Filter.u);
  //connect(PCA2.u,PCA2Filter.u);
  //connect(PCINL.u,PCINLFilter.u);

  //connect(PCA1Filter.y,w1.z1);
  //connect(PCA2Filter.y,w2.z1);
  //connect(PCINLFilter.y,p.z);

  //connect(pipelinePressureSetPoint,p.Pin_ss);

    annotation (experiment(
        StopTime=36000,
        Interval=1,
        Tolerance=1e-012,
        Algorithm="Dassl"), __Dymola_experimentSetupOutput);
  end net_PIDstruc2;

  model net_PIDstruc2_Simulate

    replaceable package components = networkModels.networkComponents;

    import SI = Modelica.SIunits;
    import Modelica.Constants.R;
    import Modelica.Constants.pi;

    networkModels.gasliftNetworks.net_PIDstruc2 netPid;

    //parameter Real[3] PresssureSetpoints = { 2.3244e+006, 2.55069e+006, 2.56983e+006};

    //2.54381000e+06;
    parameter SI.Pressure PCT1SetPoint= netPid.ZcT1*netPid.x0[4]*R*netPid.w1.par.T_r/((netPid.w1.V_r - (netPid.x0[5]-netPid.w1.par.rho_L*netPid.w1.par.L_bh*netPid.w1.par.S_bh)/netPid.w1.par.rho_L)*netPid.w1.par.M_G_r_t)
                                        - netPid.y0PCINL - netPid.m.Fric;

    //2.57117000e+06;
    parameter SI.Pressure PCT2SetPoint= netPid.ZcT2*netPid.x0[9]*R*netPid.w2.par.T_r/((netPid.w2.V_r - (netPid.x0[10]-netPid.w2.par.rho_L*netPid.w2.par.L_bh*netPid.w2.par.S_bh)/netPid.w2.par.rho_L)*netPid.w2.par.M_G_r_t)
                                          - netPid.y0PCINL - netPid.m.Fric;

    //2.27915000e+06;
    parameter SI.Pressure PCTOPSetPoint = netPid.ZcTOP*netPid.x0[14]*R*netPid.p.par.T2/(netPid.p.par.M_Gr*((pi * netPid.p.par.r2 ^ 2) * (netPid.p.par.L2 + netPid.p.par.L3) - netPid.x0[15] / netPid.p.par.rho_L)) - netPid.sep.Psep_nom;

     components.gasManifoldStep gm1(offset = 1,  height = 0.05, startTime=1*3600, fout(p(start=86.5e5)));
     components.gasManifoldStep gm2(offset = 1,  height = 0.05, startTime=1*3600, fout(p(start=91.6e5)));
     //Modelica.Blocks.Sources.Step gm1_step(offset = 1,  height = 0.07, startTime=4*3600);
     //Modelica.Blocks.Sources.Step gm2_step(offset = 1,  height = 0.07, startTime=4*3600);

     //Modelica.Blocks.Sources.Constant up(k = PCTOPSetPoint);
     Modelica.Blocks.Sources.Step step1(offset = PCTOPSetPoint,  height = -0.05e5, startTime=6*3600);
     Modelica.Blocks.Sources.Step step2(offset = 0,  height = -0.05e5, startTime=9*3600);

     //components.gasManifold gm_1;
     //components.gasManifold gm_2;

     Modelica.Blocks.Math.Add u3(u1=step1.y, k1=1, u2=step2.y, k2=1);

     Modelica.Blocks.Sources.Step u4(offset = PCT1SetPoint,  height = -0.1e5, startTime=5*3600);
     Modelica.Blocks.Sources.Step u5(offset = PCT2SetPoint,  height = -0.1e5, startTime=7*3600);

     Modelica.Blocks.Sources.Constant GOR1(k = 0);
     Modelica.Blocks.Sources.Constant GOR2(k = 0);

     //Modelica.Blocks.Sources.Constant Pr1(k = 160);
     //Modelica.Blocks.Sources.Constant Pr2(k = 170);

     Modelica.Blocks.Sources.Constant d_PI1(k = 0);
     Modelica.Blocks.Sources.Constant d_PI2(k = 0);

     Modelica.Blocks.Sources.Constant d_Psep(k = 0);

  equation
     connect(GOR1.y, netPid.w1.GOR);
     connect(GOR2.y, netPid.w2.GOR);
     connect(d_PI1.y, netPid.w1.d_PI);
     connect(d_PI2.y, netPid.w2.d_PI);
     //connect(Pr1.y, netPid.w1.P_res);
     //connect(Pr2.y, netPid.w2.P_res);
     connect(gm1.fout, netPid.w1.fin);
     connect(gm2.fout, netPid.w2.fin);
     connect(u3.y, netPid.valveDP_top_SetPoint);
     connect(u4.y, netPid.valveDP_whd1_SetPoint);
     connect(u5.y, netPid.valveDP_whd2_SetPoint);
     connect(d_Psep.y,netPid.sep.d_Psep)
    annotation (experiment(
        StopTime=36000,
        Interval=10,
        Tolerance=1e-008,
        Algorithm="Dassl"), __Dymola_experimentSetupOutput);
     connect(d_Psep.y,netPid.sep.d_Psep)
    annotation (experiment(
        StopTime=36000,
        Interval=10,
        Tolerance=1e-008,
        Algorithm="Dassl"), __Dymola_experimentSetupOutput);

    annotation (experiment(
        StopTime=36000,
        Interval=10,
        Tolerance=1e-008,
        Algorithm="Dassl"));
  end net_PIDstruc2_Simulate;
end gasliftNetworks;
