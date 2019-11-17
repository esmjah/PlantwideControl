# -*- coding: utf-8 -*-
import sys
import casadi as Ca
import control
import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as NP
from scipy import signal
import modelCompiler
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import SingleShooting
import MultipleShooting
import DLQRControl
import StateEstimator
import DaeModel
import Olga_sim as Olga
import ConfigParser

# import ConfigParser

model = 'OlgaNet.ServerDemo.'
# serverIP = 'olga.itk.ntnu.no'
serverIP = 'localhost'
NP.set_printoptions(precision=3)
reload(StateEstimator)
reload(SingleShooting)
reload(MultipleShooting)
reload(DLQRControl)
reload(Olga)
reload(DaeModel)

# measurementTagsModelica = ['net.w1.P_bh','net.w2.P_bh','net.w1.P_r_t','net.w2.P_r_t','net.w1.deltaP_valve','net.w2.deltaP_valve','net.p.deltaP_valve']
# measurementTagsOlga =     ['BTH1.PT','BTH2.PT','WHD1.PT','WHD2.PT','Prod_Choke_1.VALVDP','Prod_Choke_2.VALVDP','Valve.VALVDP']
measurementTagsModelica = ['net.w1.w_L_out', 'net.w1.w_G_out', 'net.w2.w_L_out', 'net.w2.w_G_out', 'net.w1.P_r_t',
                           'net.w2.P_r_t', 'net.w1.deltaP_valve', 'net.w2.deltaP_valve', 'net.p.deltaP_valve']
measurementTagsOlga = ['WHD1.GLTHL', 'WHD1.GG', 'WHD2.GLTHL', 'WHD2.GG', 'WHD1.PT', 'WHD2.PT', 'Prod_Choke_1.VALVDP',
                       'Prod_Choke_2.VALVDP', 'Valve.VALVDP']

# measurementTagsModelica = ['net.w1.fin.p','net.w2.fin.p','net.w1.w_L_out','net.w1.w_G_out','net.w2.w_L_out','net.w2.w_G_out','net.p.w_L_out','net.p.w_G_out','net.w1.P_r_t','net.w2.P_r_t','net.w1.deltaP_valve','net.w2.deltaP_valve','net.p.deltaP_valve']
# measurementTagsOlga =     ['ANN1.PT','ANN2.PT','WHD1.GLTHL','WHD1.GG','WHD2.GLTHL','WHD2.GG','TOP.GLTHL','TOP.GG','WHD1.PT','WHD2.PT','Prod_Choke_1.VALVDP','Prod_Choke_2.VALVDP','Valve.VALVDP']


measurementsMonitorModelica = ['net.w1.fin.p', 'net.w2.fin.p', 'net.w1.P_bh_f', 'net.w2.P_bh_f', 'net.p.w_L_out',
                               'net.p.w_G_out', 'net.PCA1.u', 'net.PCA2.u', 'net.PCINL.u', 'net.p.fin.p',
                               'net.w1.alpha_G_m_in', 'net.w2.alpha_G_m_in']
measurementsMonitorOlga = ['ANN1.PT', 'ANN2.PT', 'BTH1.PT', 'BTH2.PT', 'TOP.GLTHL', 'TOP.GG', 'PC-A1.CONTR',
                           'PC-A2.CONTR', 'PC-INL.CONTR', 'PC-INL.MEASVAR', 'IPR_1.GASFRACTION', 'IPR_2.GASFRACTION']

controlTagsModelica = ['inputGas1', 'inputGas2', 'uTOP', 'uWHD1', 'uWHD2']
controlTagsOlga = ['Gas_Source_1.MASSFLOW', 'Gas_Source_2.MASSFLOW', 'PC-TOP.SETPOINT', 'PC-T1.SETPOINT',
                   'PC-T2.SETPOINT']

openLoopControlTagsModelica = ['P_res1', 'P_res1', 'inputGas1', 'inputGas2', 'net.p.z', 'net.w1.z1', 'net.w2.z1']
openLoopControlTagsOlga = ['P_res1', 'P_res1',
                           'Gas_Source_1.MASSFLOW',
                           'Gas_Source_2.MASSFLOW',
                           'PC-INL.CONTR',
                           'PC-A1.CONTR',
                           'PC-A2.CONTR']

models_list = modelCompiler.getModelList('.', True)
ocp = modelCompiler.getOCP('ocpGaslift_GOR', models_list)

ocp.makeSemiExplicit()
ocp.eliminateIndependentParameters()
ocp.eliminateDependentParameters()
ocp.eliminateDependentParameterInterdependencies()
ocp.eliminateAlgebraic()

DT = 10  # for simulation and kalman filtering

scaleXModelica = Ca.vertcat([ocp.variable(ocp.x[k].getName()).nominal for k in range(ocp.x.size())])
scaleZModelica = Ca.vertcat([ocp.variable(ocp.z[k].getName()).nominal for k in range(ocp.z.size())])
scaleUModelica = Ca.vertcat([ocp.variable(ocp.u[k].getName()).nominal for k in range(ocp.u.size())])

## Start Kalman Filter
u = Ca.vertcat([ocp.variable(xit.getName()).initialGuess.getValue() for xit in ocp.u])

daeModel = DaeModel.DaeModel()
daeModel.init(ocp, DT, measurementTagsModelica)

x0, z0, y_0 = daeModel.findSteadyState(u, None, None, 0, consList=['net.w1.z1', 'net.w2.z1', 'net.p.z'])

uOpt, xOpt, zOpt, yOpt = daeModel.findOptimalPoint(u, x0=x0, z0=z0, simCount=0,
                                                   consList=['net.w1.z1', 'net.w2.z1', 'net.p.z'])

x0 = xOpt
z0 = zOpt
y0 = yOpt

print 'uOpt:', uOpt

measurementTagsOpt = ['net.w1.fin.p', 'net.w1.P_r_t', 'net.w1.deltaP_valve', 'net.w1.w_G_out', 'net.w1.w_L_out',
                      'GOR_hw1', 'net.w2.fin.p', 'net.w2.P_r_t', 'net.w2.deltaP_valve', 'net.w2.w_G_out',
                      'net.w2.w_L_out', 'GOR_hw2', 'net.p.fin.p', 'net.p.P2_t', 'net.p.deltaP_valve', 'net.p.w_G_out',
                      'net.p.w_L_out', 'net.p.w_mix_out', 'Obj']

measurementsEst = Ca.vertcat([ocp.variable(varName).v for varName in measurementTagsOpt])
measure_funcEst = ocp.beq(measurementsEst)
G = Ca.SXFunction(Ca.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t), [measure_funcEst])
G.init()
G.setInput(xOpt, 'x')
G.setInput(zOpt, 'z')
G.setInput(uOpt, 'p')
G.setInput(0, 't')
G.evaluate()
y_Opt = G.getOutput()
OptObj = y_Opt[-1]
print 'yOpt:', y_Opt

sys.stdout.flush()

## jacobian calculation Test

sysIn = Ca.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t)
sysOut = Ca.daeOut(ode=ocp.ode(ocp.x), alg=ocp.alg)
odeF = Ca.SXFunction(sysIn, sysOut)
odeF.init()

measTags = ['Obj']
mSX = Ca.vertcat([ocp.variable(measTags[k]).beq for k in range(len(measTags))])
mSXF = Ca.SXFunction(sysIn, [mSX])
mSXF.init()

AS = odeF.jac('x', 'ode')
BS = odeF.jac('p', 'ode')
CS = mSXF.jac('x', 0)
DS = mSXF.jac('p', 0)

funcJacs = Ca.SXFunction(sysIn, [AS, BS, CS, DS])
funcJacs.init()

scaleX = Ca.vertcat([ocp.variable(ocp.x[k].getName()).nominal for k in range(ocp.x.size())])
scaleZ = Ca.vertcat([ocp.variable(ocp.z[k].getName()).nominal for k in range(ocp.z.size())])
scaleU = Ca.vertcat([ocp.variable(ocp.u[k].getName()).nominal for k in range(ocp.u.size())])
scaleY = Ca.vertcat([ocp.variable(k).nominal for k in measurementTagsModelica])

monitorScale = Ca.vertcat(
    [ocp.variable(measurementsMonitorModelica[k]).nominal for k in range(len(measurementsMonitorModelica))])

uMin = Ca.vertcat([ocp.variable(ocp.u[k].getName()).min.getValue() for k in range(ocp.u.size())])
uMax = Ca.vertcat([ocp.variable(ocp.u[k].getName()).max.getValue() for k in range(ocp.u.size())])

# sanity check
controlTagsOCP = [ocp.u[k].getName() for k in range(ocp.u.size())]
for k in range(len(controlTagsModelica)):
    if controlTagsModelica[k] != controlTagsOCP[k]:
        raise AssertionError

qPid = 1
rPid = 1

# ocp.x =  [net.PCINL.intError, net.PCTOP.intError, net.PCA1.intError, net.PCT1.intError, net.PCA2.intError, net.PCT2.intError, net.PCT1Filter.y, net.PCT2Filter.y
# , net.PCTOPFilter.y, net.PCA1Filter.y, net.PCA2Filter.y, net.PCINLFilter.y, net.w1.m_Ga, net.w1.m_Gw, net.w1.m_Lw, net.w2.m_Ga, net.w2.m_Gw, net.w2.m_Lw, net.p.m_gp, net.p.m_lp, net.p.m_gr, net.p.m_lr]

Q = NP.diag([qPid, qPid, qPid, qPid, qPid, qPid,
             #             net.w1.m_Ga, net.w1.m_Gw, net.w1.m_Lw, net.w2.m_Ga, net.w2.m_Gw, net.w2.m_Lw, net.p.m_gp, net.p.m_lp, net.p.m_gr, net.p.m_lr]
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) * 0.1

# ['WHD1.GLTHL','WHD1.GG','WHD2.GLTHL','WHD2.GG','WHD1.PT','WHD2.PT','Prod_Choke_1.VALVDP','Prod_Choke_2.VALVDP','Valve.VALVDP']
R0 = NP.diag([100, 100, 100, 100, 1e5, 1e5, 1e4, 1e4, 1e4])  # same as NMPC
# R0 = NP.diag([1, 1, 1, 1, 1, 1, 1, 1, 1])*100 # same as old one
## Start Kalman Filter
KF = StateEstimator.StateEstimator()
KF.initEstimator(ocp, DT, measurementTagsModelica, Q=Q, R=R0)
u = Ca.vertcat([ocp.variable(xit.getName()).initialGuess.getValue() for xit in ocp.u])
x0F = xOpt  # + [NP.random.normal()*scaleX[k]/10 for k in range(ocp.x.size())]
KF.setState(x0F, z0)
# KF.computeKalmanSystem(u,x0)
# flush print buffer
print '----'
sys.stdout.flush()

# flush print buffer
print '----'
sys.stdout.flush()

## Start Simulator
openLoop = False
controlOlga = True
NIT = int(NP.ceil((3600 * 55 + 100) / DT))
xSh = []
if not openLoop:
    if controlOlga:
        da = Olga.OLGA_connect(serverIP=serverIP, modelPrefix=model)
        # Olga.OLGA_restart_time(da,True)
    else:
        xSim = NP.copy(x0) + Ca.vertcat([NP.random.normal() * scaleX[k] / 10 for k in range(ocp.x.size())])
        zSim = NP.copy(z0)
        xSh.append(xSim)

else:
    readPreviousSimulation = True
    if readPreviousSimulation:
        # dataFile ='Network-march16.mat'
        dataFile = 'yerrSim.mat'
        data = sio.loadmat(dataFile)

        yh = data['yh']
        olgaTags = data['olgaTags']

        olgaTags = [str(ti).rstrip() for ti in olgaTags]
        OlgaTagIndex = {}

        for k in range(len(olgaTags)):
            OlgaTagIndex[olgaTags[k]] = k

        NIT = min(int(NP.floor(yh[-1][OlgaTagIndex["TIME"]] / DT)), NIT)

        OlgaF = {}
        for k in olgaTags:
            OlgaF[k] = interp1d([yi[OlgaTagIndex['TIME']] for yi in yh], [yi[OlgaTagIndex[k]] for yi in yh])

    else:
        mh = []
        time = []
        x0k = x0 + [NP.random.normal() * scaleX[k] / 10 for k in range(ocp.x.size())]
        xSh.append(x0k)
        zf = None
        for k in range(NIT):
            [x0k, zf, y] = daeModel.oneStep(x0k, u, zf)
            mh.append(y)
            time.append((k + 1) * DT)
            xSh.append(x0k)

        OlgaF = {}
        for k in range(len(measurementTagsOlga)):
            OlgaF[measurementTagsOlga[k]] = interp1d(time, [yi[k] for yi in mh])

        for k in range(len(controlTagsOlga)):
            OlgaF[controlTagsOlga[k]] = interp1d([0] + time, [u[k]] + [u[k] for j in range(NIT)])

        plotSimulation = True
        if plotSimulation:
            for sit in range(len(measurementTagsOlga)):
                plt.figure()
                plt.plot([(ti + 1) * DT for ti in range(NIT)],
                         [OlgaF[measurementTagsOlga[sit]](ti + 1) * DT for ti in range(NIT)],
                         label=measurementTagsModelica[sit])
                plt.legend()

xh = []
uh = []
yh = []
yPh = []
Juh_OnSOC = []
simTime = []
monitoredOlgah = []
monitoredModelicah = []
Objh_OnSOC = []
ObjTotalh = []
xh.append(Ca.DMatrix(x0))
ULog = None

monitor = True
if monitor:
    measurementsEst = Ca.vertcat([ocp.variable(varName).v for varName in measurementsMonitorModelica])
    measure_funcEst = ocp.beq(measurementsEst)
    HEst = Ca.SXFunction(Ca.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t), [measure_funcEst])
    HEst.init()


if openLoop or controlOlga:  # testing state estimation on
    if controlOlga and not openLoop:
        y0Olga = Olga.OLGA_read_tags(da, measurementTagsOlga)

extendedKalman = True


KF.setState(x0, z0)
x_hat = KF.getState()
z_hat = Ca.vertcat([ocp.variable(ocp.z[k].getName()).start for k in range(ocp.z.size())])



# MPC = False
nX = ocp.x.size()
xMin = NP.zeros(nX)
xMax = NP.zeros(nX)

KalmanStopk = 2 * NIT

alpha_Gm1 = 0.0
alpha_Gm2 = 0.0
d_alpha = 0.03/(10.0*3600.0)
Obj_total = 0.0

gasPrice = ocp.variable('gasPrice').start


u_k = NP.array(NP.zeros(ocp.u.size()), NP.float64)
u_k_1 = NP.array(NP.zeros(ocp.u.size()), NP.float64)
u0 = NP.array(NP.zeros(ocp.u.size()), NP.float64)
NP.copyto(u_k, NP.reshape(NP.array(uOpt), ocp.u.size()))
NP.copyto(u_k_1, NP.reshape(NP.array(uOpt), ocp.u.size()))
NP.copyto(u0, NP.reshape(NP.array(uOpt), ocp.u.size()))
k0 = int(0)
k0Control = 3600 / DT


uk = Ca.DMatrix(u_k)


IntegralError = NP.array(NP.zeros(ocp.u.size()), NP.float64)
Error_k = NP.array(NP.zeros(ocp.u.size()), NP.float64)
Error_k_1 = NP.array(NP.zeros(ocp.u.size()), NP.float64)

Kp = 0.05*NP.array([1, 1, 0.2, 3, 3], NP.float64)
scale_U = NP.transpose(NP.array(scaleUModelica))[0]
Ti = 14
Td = 0
du_max = NP.array([3e-4,   3e-4,   1.0e2,   2.0e2,   2.0e2])


for k in range(k0, NIT+1):
    sys.stdout.flush()

    if (k > 3600 * 10 / DT) and (k <= 3600 * 20 / DT):
        KF.R = R0 - 0.0275 * (k - 10*3600 / DT) * NP.diag(NP.ones(len(measurementTagsOlga)))

    if (k > 3600 * 10 / DT) and (k <= 3600 * 20 / DT):
        alpha_Gm1 += d_alpha * DT
    if (k > 3600 * 30 / DT) and (k <= 3600 * 40 / DT):
        alpha_Gm2 += d_alpha * DT

    funcJacs.setInput(x_hat, 'x')
    funcJacs.setInput(uk, 'p')

    funcJacs.evaluate()

    a = funcJacs.getOutput(0)
    b = funcJacs.getOutput(1)
    c = funcJacs.getOutput(2)
    d = funcJacs.getOutput(3)

    sys_nm = control.StateSpace(a, b, c, d)
    sys_m = sys_nm.minreal()
    A = sys_m.A
    B = sys_m.B
    C = sys_m.C
    D = sys_m.D

    A_inv = NP.linalg.inv(-A)
    Ju = NP.array(C * (A_inv * B) + D)[0]

    if k >= k0Control:
        Error_k = -Ju * scale_U
        IntegralError += Error_k
        d_Error_dt = (Error_k - Error_k_1) / DT
        uk_PID = u0 + Kp * scale_U * (Error_k + IntegralError / Ti + Td * d_Error_dt)

        du = uk_PID - u_k_1

        for i in range(ocp.u.size()):
            if NP.abs(du[i]) < du_max[i]:
                u_k[i] = uk_PID[i]
            else:
                u_k[i] = u_k_1[i] + NP.sign(du[i]) * du_max[i]

            if u_k[i] < uMin[i]:
                u_k[i] = uMin[i]
            if u_k[i] > uMax[i]:
                u_k[i] = uMax[i]

        u_k_1 = u_k
        Error_k_1 = Error_k

    uk = Ca.DMatrix([u_k[0], u_k[1], u_k[2], u_k[3], u_k[4]])
    Olga.OLGA_write_tags(da, dict(zip(controlTagsOlga, uk)))

    Olga.OLGA_simulate_step(da, stepsize=DT)
    y = Olga.OLGA_read_tags(da, measurementTagsOlga)

    fOut = KF.kalman(y, uk)

    x_hat = fOut['x_hat']
    z_hat = fOut['zh']
    yP = fOut['yP']

    if monitor:
        HEst.setInput(x_hat, 'x')
        HEst.setInput(z_hat, 'z')
        HEst.setInput(uk, 'p')
        HEst.setInput((k + 1) * DT, 't')
        HEst.evaluate()

        yMonitor = Olga.OLGA_read_tags(da, measurementsMonitorOlga)

        yMonitor[-2] = alpha_Gm1
        yMonitor[-1] = alpha_Gm2

        Obj_k = gasPrice * (uk[0] + uk[1]) - (y[0] + y[2])
        Obj_total = Obj_total + Obj_k

        doPrint = True
        if doPrint and NP.mod(k * DT, 360.0) == 0:
            print "====================================================="
            print 'Simulation time:', k * DT / (3600.0), 'hours'
            print 'PredictedMeasurmentError', (y - yP)  # /KF.measurementScaling
            print 'Control Input', uk  # /scaleU
            print 'Objective', Ca.DMatrix([Obj_k, Obj_total])
            print 'Gradients', Ca.DMatrix(Ju)
            # print 'steadyStates',xOpt/scaleX
            # print '      states',x_hat
            # print 'steadyAlgStates',zOpt/scaleZ
            # print 'Alg      States',z_hat/scaleZ
            print 'Monitor    Olga', Ca.DMatrix(yMonitor)  # /monitorScale
            print 'MonitorModelica', HEst.getOutput()  # /monitorScale
            sys.stdout.flush()

    if k % 10 == 0:
        if monitor:
            monitoredOlgah.append(yMonitor)
            monitoredModelicah.append(HEst.getOutput())
            Objh_OnSOC.append(Obj_k)
            ObjTotalh.append(Obj_total)
        uh.append(uk)
        yPh.append(yP)
        xh.append(x_hat)
        yh.append(y)
        Juh_OnSOC.append(Ju)
        simTime.append(k * DT)

    # if(k==359):
    #    thebug

execfile('plotCurves_struc2.py')
execfile('SaveSimData_MESC_olga.py')
