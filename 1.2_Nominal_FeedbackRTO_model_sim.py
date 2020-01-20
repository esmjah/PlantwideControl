# -*- coding: utf-8 -*-
# Authors: Andres Codas, Esmaeil Jahanshahi
# Norwegian University of Science and Technology

# This simulation does not require the Olga Simulator. This is a test script to run the simulation without using Olga.
# Here, another instance of the simplified dynamic model is simulated and treated as the process.
# This script runs FeedbackRTO when process is at nominal conditions (no disturbances)
# Simulations starts from an initial point where the gas injection rate is 1 kg/s and all valves are 50% open.
# Optimization steers the process to the steady-state optimal point.
# Two Modelica models are used here:
# - "\modelCompiler\ocpGaslift_Sim.mop" used to simulate the gas-lift network.
# - "\modelCompiler\ocpGaslift_Nominal.mop" used to solve the optimal problem.

import sys
import os
import casadi as ca
import control
import numpy as np
import ConfigParser

_path = os.path.dirname(os.path.realpath(__file__))

os.chdir(_path)
sys.path.append(_path + '\\modelCompiler')
sys.path.append(_path + '\\StateEstimator')
sys.path.append(_path + '\\DaeModel')
sys.path.append(_path + '\\modelSimulator')

import modelCompiler
import ModelSimulator
import StateEstimator
import DaeModel

model = 'OlgaNet.ServerDemo.'
# serverIP = 'olga.itk.ntnu.no'
serverIP = 'localhost'
np.set_printoptions(precision=3)
reload(StateEstimator)
reload(DaeModel)
reload(ModelSimulator)

measurementTagsModelica = ['net.w1.w_L_out', 'net.w1.w_G_out', 'net.w2.w_L_out', 'net.w2.w_G_out', 'net.w1.P_r_t',
                           'net.w2.P_r_t', 'net.w1.deltaP_valve', 'net.w2.deltaP_valve', 'net.p.deltaP_valve']
measurementTagsSim = ['net.w1sim.w_L_out', 'net.w1sim.w_G_out', 'net.w2sim.w_L_out', 'net.w2sim.w_G_out', 'net.w1sim.P_r_t',
                           'net.w2sim.P_r_t', 'net.w1sim.deltaP_valve', 'net.w2sim.deltaP_valve', 'net.p.deltaP_valve']

measurementsMonitorModelica = ['net.w1.P_bh','net.w2.P_bh','net.p.w_L_out','net.p.w_G_out','net.PCA1.u','net.PCA2.u','net.PCINL.u','net.p.fin.p','net.w1.GOR_m_t','net.w2.GOR_m_t']
measurementsMonitorSim = ['net.w1sim.P_bh','net.w2sim.P_bh','net.p.w_L_out','net.p.w_G_out','net.PCA1.u','net.PCA2.u','net.PCINL.u','net.p.fin.p','net.w1sim.GOR_m_t','net.w2sim.GOR_m_t']

controlTagsModelica = ['inputGas1', 'inputGas2', 'uTOP','uWHD1', 'uWHD2']
controlTagsSim =  ['P_res1', 'P_res2','alpha_Gm1' , 'alpha_Gm2' , 'inputGas1', 'inputGas2', 'uTOP','uWHD1', 'uWHD2']

openLoopControlTagsModelica = ['alpha_G_m_in1','alpha_G_m_in2','inputGas1', 'inputGas2','net.p.z', 'net.w1.z1', 'net.w2.z1']
openLoopcontrolTagsSim =     ['alpha_G_m_in1','alpha_G_m_in2',
                               'Gas_Source_1.MASSFLOW',
                               'Gas_Source_2.MASSFLOW',
                               'PC-INL.CONTR',
                               'PC-A1.CONTR',
                               'PC-A2.CONTR']

tags = measurementTagsSim + measurementsMonitorSim

nTags = len(tags)
nyTags = len(measurementTagsModelica)

models_list = modelCompiler.getModelList('.', True)
ocp = modelCompiler.getOCP('ocpGaslift_Nominal', models_list)
ocpSim = modelCompiler.getOCP('ocpGaslift_Sim', models_list)

ocp.makeSemiExplicit()
ocp.eliminateIndependentParameters()
ocp.eliminateDependentParameters()
ocp.eliminateDependentParameterInterdependencies()
ocp.eliminateAlgebraic()

ocpSim.makeSemiExplicit()
ocpSim.eliminateIndependentParameters()
ocpSim.eliminateDependentParameters()
ocpSim.eliminateDependentParameterInterdependencies()
ocpSim.eliminateAlgebraic()

scaleXModelica = ca.vertcat([ocp.variable(ocp.x[k].getName()).nominal for k in range(ocp.x.size())])
scaleZModelica = ca.vertcat([ocp.variable(ocp.z[k].getName()).nominal for k in range(ocp.z.size())])
scaleUModelica = ca.vertcat([ocp.variable(ocp.u[k].getName()).nominal for k in range(ocp.u.size())])

DT = 10  # for simulation and kalman filtering

# initial values for model simulation
# u0Sim = ca.vertcat([ocpSim.variable(xit.getName()).start for xit in ocpSim.u])
u0Sim = ca.vertcat([160, 170, 0, 0, 1, 1, 75709.4, 251633, 276170])
daeModel = DaeModel.DaeModel()
daeModel.init(ocp, DT, measurementTagsModelica)

daeModelSim = DaeModel.DaeModel()
daeModelSim.init(ocpSim, DT, tags)

# Initial conditions for optimal problem
u0 = ca.vertcat([ocp.variable(xit.getName()).start for xit in ocp.u])
x0, z0, y0 = daeModel.findSteadyState(u0, None, None, 0, consList=['net.w1.z1', 'net.w2.z1', 'net.p.z'])
x0Sim, z0Sim, y_0Sim = daeModelSim.findSteadyState(u0Sim, None, None, 0,
                                                   consList=['net.w1sim.z1', 'net.w2sim.z1', 'net.p.z'])

# Finding Optimal Steady-state
uOpt, xOpt, zOpt, yOpt = daeModel.findOptimalPoint(u0, x0=x0, z0=z0, simCount=0,
                                                   consList=['net.w1.z1', 'net.w2.z1', 'net.p.z'])
uOptSim, xOptSim, zOptSim, yOptSim = daeModelSim.findOptimalPoint(u0Sim, x0=x0Sim, z0=z0Sim, simCount=0,
                                                                  consList=['net.w1sim.z1', 'net.w2sim.z1', 'net.p.z'])

print 'uOpt:', uOpt
print 'yOpt:', yOpt

sys.stdout.flush()

# jacobian calculation Test

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

uMin = ca.vertcat([ocp.variable(ocp.u[k].getName()).min.getValue() for k in range(ocp.u.size())])
uMax = ca.vertcat([ocp.variable(ocp.u[k].getName()).max.getValue() for k in range(ocp.u.size())])

monitorScale = ca.vertcat(
    [ocp.variable(measurementsMonitorModelica[k]).nominal for k in range(len(measurementsMonitorModelica))])

# sanity check
controlTagsOCP = [ocp.u[k].getName() for k in range(ocp.u.size())]
for k in range(len(controlTagsModelica)):
    if controlTagsModelica[k] != controlTagsOCP[k]:
        raise AssertionError

qPid = 1
rPid = 1

Q = np.diag([qPid, qPid, qPid, qPid, qPid, qPid,
             # net.w1.m_Ga, net.w1.m_Gw, net.w1.m_Lw, net.w2.m_Ga, net.w2.m_Gw, net.w2.m_Lw, net.p.m_gp, net.p.m_lp, net.p.m_gr, net.p.m_lr]
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) * 0.1

# ['WHD1.GLTHL','WHD1.GG','WHD2.GLTHL','WHD2.GG','WHD1.PT','WHD2.PT','Prod_Choke_1.VALVDP','Prod_Choke_2.VALVDP','Valve.VALVDP']
R0 = np.diag([1, 1, 1, 1, 1, 1, 1, 1, 1]) * 100
# R = np.diag(            [0.001,          0.001,       0.001,        0.001,     0.001,      0.1 ,   1e3,    1e3, rPid, rPid ]  )

## Start Kalman Filter
KF = StateEstimator.StateEstimator()
ModSim = ModelSimulator.ModelSimulator()
KF.initEstimator(ocp, DT, measurementTagsModelica, Q=Q, R=R0)
ModSim.initSimulator(ocpSim, DT, tags)
# u = uOpt #ca.vertcat([ocp.variable(xit.getName()).initialGuess.getValue() for xit in ocp.u])

KF.setState(x0, z0)
ModSim.setState(x0Sim, z0Sim)
# KF.computeKalmanSystem(u,x0)
# flush print buffer
print '----'
sys.stdout.flush()

## Start Simulator

MPC = True
Config = ConfigParser.ConfigParser()
Config.read('config.ini')

# flush print buffer
print '----'
sys.stdout.flush()

## Start Simulator
openLoop = False
controlOlga = False

xSh = []
xh = []
uh = []
yh = []
yPh = []
Juh_OnSOC = []
monitoredSimh = []
monitoredModelicah = []
Objh_OnSOC = []
simTime = []
ObjTotalh = []
xh.append(ca.DMatrix(x0))
ULog = None

xh.append(ca.DMatrix(x0))
ULog = None

monitor = True
if monitor:
    measurementsEst = ca.vertcat([ocp.variable(varName).v for varName in measurementsMonitorModelica])
    measure_funcEst = ocp.beq(measurementsEst)
    HEst = ca.SXFunction(ca.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t), [measure_funcEst])
    HEst.init()

extendedKalman = True

KF.setState(x0, z0)
x_hat = KF.getState()
# z_hat= ca.vertcat([ocp.variable(ocp.z[k].getName()).start for k in range(ocp.z.size())])

nX = ocp.x.size()
xMin = np.zeros(nX)
xMax = np.zeros(nX)
for k in range(nX):
    xMin[k] = ocp.variable(ocp.x[k].getName()).min.getValue()
    xMax[k] = ocp.variable(ocp.x[k].getName()).max.getValue()


Obj_total = 0.0
gasPrice = ocp.variable('gasPrice').start

u_0 = np.transpose(np.array(u0))[0]
u_k = np.transpose(np.array(u0))[0]
u_k_1 = np.transpose(np.array(u0))[0]

IntegralError = np.array(np.zeros(ocp.u.size()), np.float64)
Error_k = np.array(np.zeros(ocp.u.size()), np.float64)
Error_k_1 = np.array(np.zeros(ocp.u.size()), np.float64)

Kp = 0.005*np.array([1, 0.8, 0.2, 3, 3], np.float64)
scale_U = np.transpose(np.array(scaleUModelica))[0]
Ti = 14
Td = 0
du_max = np.array([3e-4,   3e-4,   1.0e2,   2.0e2,   2.0e2])

NIT = int(3600 * 36 / DT)
k0 = int(0)
k0Control = 3600 / DT

for k in range(k0, NIT+1):
    sys.stdout.flush()

    if (k > 3600 * 1 / DT) and (k <= 3600 * 11 / DT):
        KF.R = R0 - 0.027 * (k - 3600 / DT) * np.diag(np.ones(len(measurementTagsSim)))

    funcJacs.setInput(x_hat, 'x')
    funcJacs.setInput(u_k, 'p')

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
        Error_k = -Ju*scale_U
        IntegralError += Error_k
        d_Error_dt = (Error_k - Error_k_1) / DT
        uk_PID = u_0 + Kp * scale_U * (Error_k + IntegralError / Ti + Td * d_Error_dt)

        du = uk_PID - u_k_1

        for i in range(ocp.u.size()):
            if abs(du[i]) < du_max[i]:
                u_k[i] = uk_PID[i]
            else:
                u_k[i] = u_k_1[i] + np.sign(du[i]) * du_max[i]

            if u_k[i] < uMin[i]: u_k[i] = uMin[i]
            if u_k[i] > uMax[i]: u_k[i] = uMax[i]

        u_k_1 = u_k
        Error_k_1 = Error_k

    uk_sim = np.concatenate((np.array([160, 170, 0, 0]), u_k), axis=0)
    simOut = ModSim.Model_sim(uk_sim)
    x_k = simOut['x_k']
    y_k = simOut['y_k']
    # Olga.OLGA_simulate_step(da,stepsize=DT)
    y = y_k[0:nyTags]  # Olga.OLGA_read_tags(da,measurementTagsSim)
    yMonitor = y_k[nyTags:nTags]  # Olga.OLGA_read_tags(da,measurementsMonitorSim)
    Obj_k = gasPrice * (u_k[0] + u_k[1]) - (y[0] + y[2])
    Obj_total = Obj_total + Obj_k

    fOut = KF.kalman(np.array(y), u_k)
    x_hat = fOut['x_hat']
    z_hat = fOut['zh']
    yP = fOut['yC']

    HEst.setInput(x_hat, 'x')
    HEst.setInput(z_hat, 'z')
    HEst.setInput(u_k, 'p')
    HEst.setInput((k + 1) * DT, 't')
    HEst.evaluate()

    doPrint = True
    if doPrint and k % 36 == 0:  # and np.mod((k + 1 - k0MPC) * DT, DTMPC / 10.0) == 0:
        print "====================================================="
        print 'Simulation time:', k * DT / 3600.0, 'hours'
        print 'PredictedMeasurmentError', (y - yP)  # /KF.measurementScaling
        print 'OptimControl', u_k  # /scaleU
        print 'Objective', ca.DMatrix([Obj_k, Obj_total])
        print 'Gradients', ca.DMatrix(Ju)

        # print 'steadyStates',xOpt/scaleX
        print '      states', x_hat
        # print 'steadyAlgStates',zOpt/scaleZ
        # print 'Alg      States',z_hat/scaleZ
        print 'Monitor    Olga', (yMonitor)  # /monitorScale
        print 'MonitorModelica', HEst.getOutput()  # /monitorScale
        sys.stdout.flush()

    if k % 10 == 0:
        if monitor:
            monitoredSimh.append(yMonitor)
            monitoredModelicah.append(HEst.getOutput())
            Objh_OnSOC.append(Obj_k)
            ObjTotalh.append(Obj_total)
        uh.append(ca.DMatrix(u_k.tolist()))
        yPh.append(yP)
        xh.append(x_hat)
        yh.append(y)
        Juh_OnSOC.append(Ju)
        simTime.append(k * DT)

execfile('SavedResults\\plotCurves_model.py')