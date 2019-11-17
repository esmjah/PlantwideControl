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

measurementsMonitorModelica = ['net.w1.fin.p', 'net.w2.fin.p', 'net.w1.P_bh', 'net.w2.P_bh', 'net.p.w_L_out',
                               'net.p.w_G_out', 'net.PCA1.u', 'net.PCA2.u', 'net.PCINL.u', 'net.p.fin.p',
                               'net.w1.w_G_out', 'net.w2.w_G_out']
measurementsMonitorOlga = ['ANN1.PT', 'ANN2.PT', 'BTH1.PT', 'BTH2.PT', 'TOP.GLTHL', 'TOP.GG', 'PC-A1.CONTR',
                           'PC-A2.CONTR', 'PC-INL.CONTR', 'PC-INL.MEASVAR', 'WHD1.GG','WHD2.GG']

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
ocp = modelCompiler.getOCP('ocpGasliftNet', models_list)

ocp.makeSemiExplicit()
ocp.eliminateIndependentParameters()
ocp.eliminateDependentParameters()
ocp.eliminateDependentParameterInterdependencies()
ocp.eliminateAlgebraic()

DT = 10  # for simulation and kalman filtering
DTMPC = 3600  ## for the MPC algorithm  Please always choose one to be multiple of the other

scaleXModelica = Ca.vertcat([ocp.variable(ocp.x[k].getName()).nominal for k in range(ocp.x.size())])
scaleZModelica = Ca.vertcat([ocp.variable(ocp.z[k].getName()).nominal for k in range(ocp.z.size())])
scaleUModelica = Ca.vertcat([ocp.variable(ocp.u[k].getName()).nominal for k in range(ocp.u.size())])

## Start Kalman Filter
u = Ca.vertcat([ocp.variable(xit.getName()).initialGuess.getValue() for xit in ocp.u])

daeModel = DaeModel.DaeModel()
daeModel.init(ocp, DT, measurementTagsModelica)

x0, z0, y0 = daeModel.findSteadyState(u, None, None, 0, consList=['net.w1.z1', 'net.w2.z1', 'net.p.z'])

uOpt, xOpt, zOpt, yOpt = daeModel.findOptimalPoint(u, x0=x0, z0=z0, simCount=0,
                                                   consList=['net.w1.z1', 'net.w2.z1', 'net.p.z'])

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
             # net.w1.m_Ga, net.w1.m_Gw, net.w1.m_Lw, net.w2.m_Ga, net.w2.m_Gw, net.w2.m_Lw, net.p.m_gp, net.p.m_lp, net.p.m_gr, net.p.m_lr]
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) * 0.1

# ['WHD1.GLTHL','WHD1.GG','WHD2.GLTHL','WHD2.GG','WHD1.PT','WHD2.PT','Prod_Choke_1.VALVDP','Prod_Choke_2.VALVDP','Valve.VALVDP']
R0 = NP.diag([1, 1, 1, 1, 1, 1, 1, 1, 1]) * 100
# Start Kalman Filter
KF = StateEstimator.StateEstimator()
KF.initEstimator(ocp, DT, measurementTagsModelica, Q=Q, R=R0)
u = Ca.vertcat([ocp.variable(xit.getName()).initialGuess.getValue() for xit in ocp.u])
x0F = xOpt  # + [NP.random.normal()*scaleX[k]/10 for k in range(ocp.x.size())]
KF.setState(x0F, z0)
# KF.computeKalmanSystem(u,x0)

# flush print buffer
print '----'
sys.stdout.flush()

# Start Simulator
openLoop = False
controlOlga = True
if controlOlga:
    da = Olga.OLGA_connect(serverIP=serverIP, modelPrefix=model)

xSh = []
xh = []
uh = []
yh = []
simTime = []
yPh = []
Juh_SOC = []
monitoredOlgah = []
monitoredModelicah = []
Objh = []
ObjTotalh = []
y1h = []
y2h = []
xh.append(Ca.DMatrix(x0))
ULog = None

monitor = True
if monitor:
    measurementsEst = Ca.vertcat([ocp.variable(varName).v for varName in measurementsMonitorModelica])
    measure_funcEst = ocp.beq(measurementsEst)
    HEst = Ca.SXFunction(Ca.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t), [measure_funcEst])
    HEst.init()

    A_g = ocp.beq('net.p.A_g')
    A_l = ocp.beq('net.p.A_l')
    A_1 = ocp.beq('net.p.A1')

    h1 = ocp.beq('net.p.h1')
    h1ss = ocp.beq('net.p.h1ss')
    hc = ocp.beq('net.p.hc')

    Alpha_Lt = ocp.beq('net.p.Alpha_Lt')
    Alpha_L2_av = ocp.beq('net.p.Alpha_L2_av')

    check = Ca.SXFunction(Ca.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t), [
        Ca.vertcat([A_g / A_1 * 100, A_l / A_1 * 100, h1 / hc * 100, h1ss / hc * 100, Alpha_Lt, Alpha_L2_av])])
    check.init()

if openLoop or controlOlga:  # testing state estimation on
    if controlOlga and not openLoop:
        y0Olga = Olga.OLGA_read_tags(da, measurementTagsOlga)

extendedKalman = True
# u0 =  Ca.vertcat([ocp.variable(xit.getName()).initialGuess.getValue() for xit in ocp.u])
# x0,z0,y0 = daeModel.findSteadyState(u0,None,None,0,consList=['net.w1.z1','net.w2.z1','net.p.z'])


KF.setState(x0, z0)
x_hat = KF.getState()
z_hat = Ca.vertcat([ocp.variable(ocp.z[k].getName()).start for k in range(ocp.z.size())])

# MPC = False
nX = ocp.x.size()
xMin = NP.zeros(nX)
xMax = NP.zeros(nX)
for k in range(nX):
    xMin[k] = ocp.variable(ocp.x[k].getName()).min.getValue()
    xMax[k] = ocp.variable(ocp.x[k].getName()).max.getValue()

NIT = int(NP.ceil(3600 * 36 / DT))
KalmanStopk = 2 * NIT

Obj_total = 0.0
gasPrice = ocp.variable('gasPrice').start


u_k = NP.array(NP.zeros(ocp.u.size()), NP.float64)
u_k_1 = NP.array(NP.zeros(ocp.u.size()), NP.float64)
u0 = NP.array(NP.zeros(ocp.u.size()), NP.float64)
NP.copyto(u_k, NP.reshape(NP.array(u), ocp.u.size()))
NP.copyto(u_k_1, NP.reshape(NP.array(u), ocp.u.size()))
NP.copyto(u0, NP.reshape(NP.array(u), ocp.u.size()))
uk = u
k0 = int(0)
k0Control = 3600 / DT


IntegralError = NP.array(NP.zeros(ocp.u.size()), NP.float64)
Error_k = NP.array(NP.zeros(ocp.u.size()), NP.float64)
Error_k_1 = NP.array(NP.zeros(ocp.u.size()), NP.float64)

Kp = NP.array([0.025, 0.025, 5e-7, 5e-7, 5e-7], NP.float64)
scale_U = NP.transpose(NP.array(scaleUModelica))[0]
Ti = 300
Td = 0
du_max = NP.array([3e-4,   3e-4,   1.0e2,   2.0e2,   2.0e2])

r_WG1 = 1.305436
r_WG2 = 1.3344603


y_WG1 = r_WG1
y_WG2 = r_WG2


for k in range(k0, NIT):
    sys.stdout.flush()

    if (k > 3600 * 1 / DT) and (k <= 3600 * 11 / DT):
        KF.R = R0 - 0.027 * (k - 3600 / DT) * NP.diag(NP.ones(len(measurementTagsOlga)))

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
        Error_k = NP.array([r_WG1-y_WG1, r_WG2-y_WG2, uMin[2]-uk[2], uMin[3]-uk[3], uMin[4]-uk[4]])
        IntegralError += Error_k*DT
        d_Error_dt = (Error_k - Error_k_1) / DT
        uk_PID = u0 + Kp * scale_U * (Error_k + IntegralError / Ti + Td * d_Error_dt)

        du = uk_PID - u_k_1

        for i in range(ocp.u.size()):
            if NP.abs(du[i]) < du_max[i]:
                u_k[i] = uk_PID[i]
            else:
                u_k[i] = u_k_1[i] + NP.sign(du[i]) * du_max[i]

            if u_k[i] < uMin[i]: u_k[i] = uMin[i]
            if u_k[i] > uMax[i]: u_k[i] = uMax[i]

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

        check.setInput(x_hat, 'x')
        check.setInput(z_hat, 'z')
        check.setInput(uk, 'p')
        check.setInput((k + 1) * DT, 't')
        check.evaluate()

        yMonitor = Olga.OLGA_read_tags(da, measurementsMonitorOlga)
        y_WG1 = yMonitor[-2]
        y_WG2 = yMonitor[-1]

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
            Objh.append(Obj_k)
            ObjTotalh.append(Obj_total)
        uh.append(uk)
        yPh.append(yP)
        xh.append(x_hat)
        yh.append(y)
        Juh_SOC.append(Ju)
        y1h.append(y_WG1)
        y2h.append(y_WG2)
        simTime.append(k * DT)

        # if(k==359):
        #    thebug

execfile('SaveSimData_SOC_olga.py')
execfile('plotCurves_struc2.py')

