# -*- coding: utf-8 -*-
# Authors: Andres Codas, Esmaeil Jahanshahi
# Norwegian University of Science and Technology

# This simulation requires the Olga Simulator. The Olga OPC Server must be running as Administrator user.
# This script runs Self-optimizing control when process is at nominal conditions (no disturbances)
# Simulations starts from an initial point where the gas injection rate is 1 kg/s and all valves are 50% open.
# Self-optimizing control steers the process to the steady-state optimal point.
# The Controlled Variavles (CV1) for Self-optimizing control are gas flow rates at the well-head ('WHD1.GG' and 'WHD2.GG').
# In practice, simplified dynamic model is not required for Self-optimizing control.
# ocpGaslift_Nominal used to estimate the gradients for comparison to other optimization methods.

import sys
import os
import casadi as ca
import control
import numpy as np

_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(_path)
sys.path.append(_path+'\\NMPC')
sys.path.append(_path+'\\modelCompiler')
sys.path.append(_path+'\\StateEstimator')
sys.path.append(_path+'\\DaeModel')
sys.path.append(_path+'\\Olga_sim')

import SingleShooting
import MultipleShooting
import modelCompiler
import StateEstimator
import DaeModel
import Olga_sim as Olga

model = 'OlgaNet.ServerDemo.'
# serverIP = 'olga.itk.ntnu.no'
serverIP = 'localhost'
np.set_printoptions(precision=3)
reload(StateEstimator)
reload(SingleShooting)
reload(MultipleShooting)
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
ocp = modelCompiler.getOCP('ocpGaslift_Nominal', models_list)

ocp.makeSemiExplicit()
ocp.eliminateIndependentParameters()
ocp.eliminateDependentParameters()
ocp.eliminateDependentParameterInterdependencies()
ocp.eliminateAlgebraic()

DT = 10  # for simulation and kalman filtering
DTMPC = 3600  ## for the MPC algorithm  Please always choose one to be multiple of the other

scaleXModelica = ca.vertcat([ocp.variable(ocp.x[k].getName()).nominal for k in range(ocp.x.size())])
scaleZModelica = ca.vertcat([ocp.variable(ocp.z[k].getName()).nominal for k in range(ocp.z.size())])
scaleUModelica = ca.vertcat([ocp.variable(ocp.u[k].getName()).nominal for k in range(ocp.u.size())])

## Start Kalman Filter
u = ca.vertcat([ocp.variable(xit.getName()).initialGuess.getValue() for xit in ocp.u])

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

measurementsEst = ca.vertcat([ocp.variable(varName).v for varName in measurementTagsOpt])
measure_funcEst = ocp.beq(measurementsEst)
G = ca.SXFunction(ca.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t), [measure_funcEst])
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

sysIn = ca.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t)
sysOut = ca.daeOut(ode=ocp.ode(ocp.x), alg=ocp.alg)
odeF = ca.SXFunction(sysIn, sysOut)
odeF.init()

measTags = ['Obj']
mSX = ca.vertcat([ocp.variable(measTags[k]).beq for k in range(len(measTags))])
mSXF = ca.SXFunction(sysIn, [mSX])
mSXF.init()

AS = odeF.jac('x', 'ode')
BS = odeF.jac('p', 'ode')
CS = mSXF.jac('x', 0)
DS = mSXF.jac('p', 0)

funcJacs = ca.SXFunction(sysIn, [AS, BS, CS, DS])
funcJacs.init()

scaleX = ca.vertcat([ocp.variable(ocp.x[k].getName()).nominal for k in range(ocp.x.size())])
scaleZ = ca.vertcat([ocp.variable(ocp.z[k].getName()).nominal for k in range(ocp.z.size())])
scaleU = ca.vertcat([ocp.variable(ocp.u[k].getName()).nominal for k in range(ocp.u.size())])
scaleY = ca.vertcat([ocp.variable(k).nominal for k in measurementTagsModelica])

monitorScale = ca.vertcat(
    [ocp.variable(measurementsMonitorModelica[k]).nominal for k in range(len(measurementsMonitorModelica))])

uMin = ca.vertcat([ocp.variable(ocp.u[k].getName()).min.getValue() for k in range(ocp.u.size())])
uMax = ca.vertcat([ocp.variable(ocp.u[k].getName()).max.getValue() for k in range(ocp.u.size())])


# sanity check
controlTagsOCP = [ocp.u[k].getName() for k in range(ocp.u.size())]
for k in range(len(controlTagsModelica)):
    if controlTagsModelica[k] != controlTagsOCP[k]:
        raise AssertionError

qPid = 1
rPid = 1

# ocp.x =  [net.PCINL.intError, net.PCTOP.intError, net.PCA1.intError, net.PCT1.intError, net.PCA2.intError, net.PCT2.intError, net.PCT1Filter.y, net.PCT2Filter.y
# , net.PCTOPFilter.y, net.PCA1Filter.y, net.PCA2Filter.y, net.PCINLFilter.y, net.w1.m_Ga, net.w1.m_Gw, net.w1.m_Lw, net.w2.m_Ga, net.w2.m_Gw, net.w2.m_Lw, net.p.m_gp, net.p.m_lp, net.p.m_gr, net.p.m_lr]

Q = np.diag([qPid, qPid, qPid, qPid, qPid, qPid,
             # net.w1.m_Ga, net.w1.m_Gw, net.w1.m_Lw, net.w2.m_Ga, net.w2.m_Gw, net.w2.m_Lw, net.p.m_gp, net.p.m_lp, net.p.m_gr, net.p.m_lr]
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) * 0.1

# ['WHD1.GLTHL','WHD1.GG','WHD2.GLTHL','WHD2.GG','WHD1.PT','WHD2.PT','Prod_Choke_1.VALVDP','Prod_Choke_2.VALVDP','Valve.VALVDP']
R0 = np.diag([1, 1, 1, 1, 1, 1, 1, 1, 1]) * 100
# Start Kalman Filter
KF = StateEstimator.StateEstimator()
KF.initEstimator(ocp, DT, measurementTagsModelica, Q=Q, R=R0)
u = ca.vertcat([ocp.variable(xit.getName()).initialGuess.getValue() for xit in ocp.u])
x0F = xOpt  # + [np.random.normal()*scaleX[k]/10 for k in range(ocp.x.size())]
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
xh.append(ca.DMatrix(x0))
ULog = None

monitor = True
if monitor:
    measurementsEst = ca.vertcat([ocp.variable(varName).v for varName in measurementsMonitorModelica])
    measure_funcEst = ocp.beq(measurementsEst)
    HEst = ca.SXFunction(ca.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t), [measure_funcEst])
    HEst.init()

    A_g = ocp.beq('net.p.A_g')
    A_l = ocp.beq('net.p.A_l')
    A_1 = ocp.beq('net.p.A1')

    h1 = ocp.beq('net.p.h1')
    h1ss = ocp.beq('net.p.h1ss')
    hc = ocp.beq('net.p.hc')

    Alpha_Lt = ocp.beq('net.p.Alpha_Lt')
    Alpha_L2_av = ocp.beq('net.p.Alpha_L2_av')

    check = ca.SXFunction(ca.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t), [
        ca.vertcat([A_g / A_1 * 100, A_l / A_1 * 100, h1 / hc * 100, h1ss / hc * 100, Alpha_Lt, Alpha_L2_av])])
    check.init()

if openLoop or controlOlga:  # testing state estimation on
    if controlOlga and not openLoop:
        y0Olga = Olga.OLGA_read_tags(da, measurementTagsOlga)

extendedKalman = True
# u0 =  ca.vertcat([ocp.variable(xit.getName()).initialGuess.getValue() for xit in ocp.u])
# x0,z0,y0 = daeModel.findSteadyState(u0,None,None,0,consList=['net.w1.z1','net.w2.z1','net.p.z'])


KF.setState(x0, z0)
x_hat = KF.getState()
z_hat = ca.vertcat([ocp.variable(ocp.z[k].getName()).start for k in range(ocp.z.size())])

# MPC = False
nX = ocp.x.size()
xMin = np.zeros(nX)
xMax = np.zeros(nX)
for k in range(nX):
    xMin[k] = ocp.variable(ocp.x[k].getName()).min.getValue()
    xMax[k] = ocp.variable(ocp.x[k].getName()).max.getValue()

NIT = int(np.ceil(3600 * 22 / DT))
KalmanStopk = 2 * NIT

Obj_total = 0.0
gasPrice = ocp.variable('gasPrice').start


u_k = np.array(np.zeros(ocp.u.size()), np.float64)
u_k_1 = np.array(np.zeros(ocp.u.size()), np.float64)
u0 = np.array(np.zeros(ocp.u.size()), np.float64)
np.copyto(u_k, np.reshape(np.array(u), ocp.u.size()))
np.copyto(u_k_1, np.reshape(np.array(u), ocp.u.size()))
np.copyto(u0, np.reshape(np.array(u), ocp.u.size()))
uk = u
u_PID_k_1 = np.transpose(np.array(u, np.float64))[0]
k0 = int(0)
k0Control = 3600 / DT


IntegralError = np.array(np.zeros(ocp.u.size()), np.float64)
Error_k = np.array(np.zeros(ocp.u.size()), np.float64)
Error_k_1 = np.array(np.zeros(ocp.u.size()), np.float64)

Kp = np.array([0.025, 0.025, 5e-7, 5e-7, 5e-7], np.float64)
scale_U = np.transpose(np.array(scaleUModelica))[0]
Ti = 150
Td = 0

du_max = np.array([6e-4,   6e-4,   4.0e2,   6.0e2,   6.0e2])

r_WG1 = 1.305436
r_WG2 = 1.3344603


y_WG1 = r_WG1
y_WG2 = r_WG2


for k in range(k0, NIT):
    sys.stdout.flush()

    if (k > 3600 * 1 / DT) and (k <= 3600 * 11 / DT):
        KF.R = R0 - 0.027 * (k - 3600 / DT) * np.diag(np.ones(len(measurementTagsOlga)))

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

    A_inv = np.linalg.inv(-A)
    Ju = np.array(C * (A_inv * B) + D)[0]

    if k >= k0Control:
        Error_k = np.array([r_WG1-y_WG1, r_WG2-y_WG2, uMin[2]-uk[2], uMin[3]-uk[3], uMin[4]-uk[4]])
        IntegralError += Error_k*DT
        d_Error_dt = (Error_k - Error_k_1) / DT
        uk_PID = u0 + Kp * scale_U * (Error_k + IntegralError / Ti + Td * d_Error_dt)

        du = uk_PID - u_k_1

        for i in range(ocp.u.size()):
            if np.abs(du[i]) <= du_max[i]:
                u_k[i] = uk_PID[i]
            else:
                u_k[i] = u_k_1[i] + np.sign(du[i]) * du_max[i]

            if u_k[i] < uMin[i]: u_k[i] = uMin[i]
            if u_k[i] > uMax[i]: u_k[i] = uMax[i]

        u_k_1 = u_k
        Error_k_1 = Error_k

    uk = ca.DMatrix([u_k[0], u_k[1], u_k[2], u_k[3], u_k[4]])
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
        if doPrint and np.mod(k * DT, 360.0) == 0:
            print "====================================================="
            print 'Simulation time:', k * DT / (3600.0), 'hours'
            print 'PredictedMeasurmentError', (y - yP)  # /KF.measurementScaling
            print 'Control Input', uk  # /scaleU
            print 'Objective', ca.DMatrix([Obj_k, Obj_total])
            print 'Gradients', ca.DMatrix(Ju)
            # print 'steadyStates',xOpt/scaleX
            # print '      states',x_hat
            # print 'steadyAlgStates',zOpt/scaleZ
            # print 'Alg      States',z_hat/scaleZ
            print 'Monitor    Olga', ca.DMatrix(yMonitor)  # /monitorScale
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

execfile('SavedResults\\SaveSimData_SelfOptimizing_olga.py')
execfile('SavedResults\\plotCurves_olga.py')
