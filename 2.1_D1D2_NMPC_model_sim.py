# -*- coding: utf-8 -*-
# Authors: Andres Codas, Esmaeil Jahanshahi
# Norwegian University of Science and Technology

# This is a test script to run the simulation without using the Olga simulator.
# This script runs NMPC in presence of disturbances in reservoir pressures (D1 and D2)
# Two Modelica models are used here:
# - ocpGaslift_Sim.mop used to simulate the gas-lift network
# - ocpGaslift_D1D2.mop used to solve the optimal problem. It has 2 augmented state to estimate the reservoir pressures

import sys
import os
import casadi as ca
import control
import numpy as np
import ConfigParser
import time

tic = time.time()

_path = os.path.dirname(os.path.realpath(__file__))

os.chdir(_path)
sys.path.append(_path + '\\NMPC')
sys.path.append(_path + '\\modelCompiler')
sys.path.append(_path + '\\StateEstimator')
sys.path.append(_path + '\\DaeModel')
sys.path.append(_path + '\\modelSimulator')

import modelCompiler
import ModelSimulator
import SingleShooting
import MultipleShooting
import StateEstimator
import DaeModel

model = 'OlgaNet.ServerDemo.'
# serverIP = 'olga.itk.ntnu.no'
serverIP = 'localhost'
np.set_printoptions(precision=3)
reload(StateEstimator)
reload(SingleShooting)
reload(MultipleShooting)
reload(DaeModel)
reload(ModelSimulator)

measurementTagsModelica = ['net.w1.w_L_out', 'net.w1.w_G_out', 'net.w2.w_L_out', 'net.w2.w_G_out', 'net.w1.P_r_t',
                           'net.w2.P_r_t', 'net.w1.deltaP_valve', 'net.w2.deltaP_valve', 'net.p.deltaP_valve']
measurementTagsSim = ['net.w1sim.w_L_out', 'net.w1sim.w_G_out', 'net.w2sim.w_L_out', 'net.w2sim.w_G_out', 'net.w1sim.P_r_t',
                           'net.w2sim.P_r_t', 'net.w1sim.deltaP_valve', 'net.w2sim.deltaP_valve', 'net.p.deltaP_valve']

measurementsMonitorModelica = ['net.w1.P_bh','net.w2.P_bh','net.p.w_L_out','net.p.w_G_out','net.PCA1.u','net.PCA2.u','net.PCINL.u','net.p.fin.p','net.w1.GOR_m_t','net.w2.GOR_m_t','net.w1.P_res','net.w2.P_res']
measurementsMonitorSim = ['net.w1sim.P_bh','net.w2sim.P_bh','net.p.w_L_out','net.p.w_G_out','net.PCA1.u','net.PCA2.u','net.PCINL.u','net.p.fin.p','net.w1sim.GOR_m_t','net.w2sim.GOR_m_t','net.w1sim.P_res','net.w2sim.P_res']

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
ocp = modelCompiler.getOCP('ocpGaslift_D1D2', models_list)
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

DT = 10  # for simulation and kalman filtering
DTMPC = 600  ## for the MPC algorithm  Please always choose one to be multiple of the other
prediction_horizon = 8 * 3600

scaleXModelica = ca.vertcat([ocp.variable(ocp.x[k].getName()).nominal for k in range(ocp.x.size())])
scaleZModelica = ca.vertcat([ocp.variable(ocp.z[k].getName()).nominal for k in range(ocp.z.size())])
scaleUModelica = ca.vertcat([ocp.variable(ocp.u[k].getName()).nominal for k in range(ocp.u.size())])

## Start Kalman Filter
u0 = ca.vertcat([ocp.variable(xit.getName()).start for xit in ocp.u])
u0Sim = ca.vertcat([ocpSim.variable(xit.getName()).start for xit in ocpSim.u])

daeModel = DaeModel.DaeModel()
daeModel.init(ocp, DT, measurementTagsModelica)

daeModelSim = DaeModel.DaeModel()
daeModelSim.init(ocpSim, DT, tags)

x0, z0, y0 = daeModel.findSteadyState(u0, None, None, 0, consList=['net.w1.z1', 'net.w2.z1', 'net.p.z'])
x0Sim, z0Sim, y_0Sim = daeModelSim.findSteadyState(u0Sim, None, None, 0,
                                                   consList=['net.w1sim.z1', 'net.w2sim.z1', 'net.p.z'])

uOpt, xOpt, zOpt, yOpt = daeModel.findOptimalPoint(u0, x0=x0, z0=z0, simCount=0,
                                                   consList=['net.w1.z1', 'net.w2.z1', 'net.p.z'])
uOptSim, xOptSim, zOptSim, yOptSim = daeModelSim.findOptimalPoint(u0Sim, x0=x0Sim, z0=z0Sim, simCount=0,
                                                                  consList=['net.w1sim.z1', 'net.w2sim.z1', 'net.p.z'])

print 'uOpt:', uOpt
print 'yOpt:', yOpt

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

# sanity check
controlTagsOCP = [ocp.u[k].getName() for k in range(ocp.u.size())]
for k in range(len(controlTagsModelica)):
    if controlTagsModelica[k] != controlTagsOCP[k]:
        raise AssertionError

qPid = 1
rPid = 1

Q = np.diag([qPid, qPid, qPid, qPid, qPid, qPid,
             #             net.w1.m_Ga, net.w1.m_Gw, net.w1.m_Lw, net.w2.m_Ga, net.w2.m_Gw, net.w2.m_Lw, net.p.m_gp, net.p.m_lp, net.p.m_gr, net.p.m_lr]
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) * 0.1

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

if MPC:
    #    controlCost = 10000*np.diag([1,1,1,1,1])
    #    finalStateCost = 0
    #    ss = SingleShooting.SingleShooting()  ## the steady state optimizer is implemented here.
    #    ss.defineOCP(ocp,DT=DTMPC,controlCost=controlCost,xOpt=xOpt,uOpt=uOpt,finalStateCost=finalStateCost)
    if Config.getboolean('MultipleShooting', 'MS'):
        solver = MultipleShooting.MultipleShooting()
    else:
        solver = SingleShooting.SingleShooting()

    deltaUCons = np.array([0.03, 0.03, 2.0e3, 5.0e3, 5.0e3]) / DTMPC  ## given in SI

    ## check consistency
    if deltaUCons.size != 0:
        assert deltaUCons.size == ocp.u.size()
    controlCost = np.diag([20.0, 20.0, 1.0, 1.0, 1.0])
    finalStateCost = 0.0
    solver.defineOCP(ocp, DT=DTMPC, controlCost=controlCost, xOpt=x0, uOpt=u0, finalStateCost=finalStateCost,
                     deltaUCons=deltaUCons)

# flush print buffer
print '----'
sys.stdout.flush()

## Start Simulator
openLoop = False
controlOlga = False

xh = []
uh = []
yh = []
yPh = []
simTime = []
Juh_NMPC = []
monitoredSimh = []
monitoredModelicah = []
Objh_NMPC = []
ObjTotalh = []

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

uk = u0
nX = ocp.x.size()
xMin = ca.DMatrix(np.zeros(nX))
xMax = ca.DMatrix(np.zeros(nX))
for ix in range(nX):
    xMin[ix] = ocp.variable(ocp.x[ix].getName()).min.getValue()
    xMax[ix] = ocp.variable(ocp.x[ix].getName()).max.getValue()

if MPC:
    ocp.variable('net.w1.P_res').min = x_hat[6]
    ocp.variable('net.w1.P_res').max = x_hat[6]
    ocp.variable('net.w2.P_res').min = x_hat[12]
    ocp.variable('net.w2.P_res').max = x_hat[12]
    for ix in [6, 12]:
        xMin[ix] = ocp.variable(ocp.x[ix].getName()).min.getValue()
        xMax[ix] = ocp.variable(ocp.x[ix].getName()).max.getValue()
    print 'xMax:', xMax[6], ',', xMax[12]
    print 'x0:', x_hat[6], ',', x_hat[12]
    print 'xMin:', xMin[6], ',', xMin[12]
    uk, objValue, stats = solver.solveProblem(x0=x_hat, u0=uk)
    print 'controlInput', uk  # /scaleU
    print 'Obj. Value:', objValue / prediction_horizon
    print 'Iterations:', stats['iter_count']
    print 'Comp. Time:', stats['t_mainloop']
    if 'return_status' in stats:
        exit_msg = stats['return_status']
    else:
        exit_msg = 'Max_CPU_time_exceeded'
    print 'Exit:', exit_msg

    if exit_msg == 'Infeasible_Problem_Detected' or exit_msg == 'Invalid_Number_Detected':
        # uk = uk_1
        print 'shit happens!'
        solver.plotNLP(x_hat, DTplot=10, additionalplottags=None)
        solver.plotFigures()
        solver.stop_me_please()
    else:
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

P_res1 = 160.0
P_res2 = 170.0
Obj_total = 0.0
dP_res = 5.0 / (10 * 3600)
gasPrice = ocp.variable('gasPrice').start

NIT = int(np.ceil(3600 * 55 / DT))
k0 = int(0)
k0MPC = 2*3600 / DT

for k in range(k0, NIT+1):
    sys.stdout.flush()

    if (k > 3600 * 10 / DT) and (k <= 3600 * 20 / DT):
        KF.R = R0 - 0.0275 * (k - 10*3600 / DT) * np.diag(np.ones(len(measurementTagsSim)))

    if (k > 3600 * 10 / DT) and (k <= 3600 * 20 / DT):
        P_res1 -= dP_res * DT
    if (k > 3600 * 30 / DT) and (k <= 3600 * 40 / DT):
        P_res2 -= dP_res * DT

    x_hat[6] = max(x_hat[6], 0)
    x_hat[12] = max(x_hat[12], 0)

    if (np.mod((k - k0MPC) * DT, DTMPC) == 0) and (k >= k0MPC):
        ocp.variable('net.w1.P_res').min = x_hat[6]
        ocp.variable('net.w1.P_res').max = x_hat[6]
        ocp.variable('net.w2.P_res').min = x_hat[12]
        ocp.variable('net.w2.P_res').max = x_hat[12]
        solver.moveHorizion()
        for ix in [6, 12]:
            xMin[ix] = ocp.variable(ocp.x[ix].getName()).min.getValue()
            xMax[ix] = ocp.variable(ocp.x[ix].getName()).max.getValue()

        print 'xMax:', xMax[6], ',', xMax[12]
        print 'x0:', x_hat[6], ',', x_hat[12]
        print 'xMin:', xMin[6], ',', xMin[12]
        uk, objValue, stats = solver.solveProblem(x_hat, uk)
        print 'controlInput', uk  # /scaleU
        print 'Obj. Value:', objValue / prediction_horizon
        print 'Iterations:', stats['iter_count']
        print 'Comp. Time:', stats['t_mainloop']
        if 'return_status' in stats:
            exit_msg = stats['return_status']
        else:
            exit_msg = 'Max_CPU_time_exceeded'
        print 'Exit:', exit_msg

        if exit_msg == 'Infeasible_Problem_Detected' or exit_msg == 'Invalid_Number_Detected':
            # uk = uk_1
            print 'shit happens!'
            solver.plotNLP(x_hat, DTplot=10, additionalplottags=None)
            solver.plotFigures()
            solver.stop_me_please()
        else:
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
    uk_sim = np.concatenate((np.array([P_res1, P_res2, 0, 0]), np.transpose(np.array(uk))[0]), axis=0)
    simOut = ModSim.Model_sim(uk_sim)
    x_k = simOut['x_k']
    y_k = simOut['y_k']
    # Olga.OLGA_simulate_step(da,stepsize=DT)
    y = y_k[0:nyTags]  # Olga.OLGA_read_tags(da,measurementTagsSim)
    yMonitor = y_k[nyTags:nTags]  # Olga.OLGA_read_tags(da,measurementsMonitorSim)
    Obj_k = gasPrice * (uk[0] + uk[1]) - (y[0] + y[2])
    Obj_total = Obj_total + Obj_k

    fOut = KF.kalman(np.array(y), uk)
    x_hat = fOut['x_hat']
    z_hat = fOut['zh']
    yP = fOut['yC']

    HEst.setInput(x_hat, 'x')
    HEst.setInput(z_hat, 'z')
    HEst.setInput(uk, 'p')
    HEst.setInput((k + 1) * DT, 't')
    HEst.evaluate()

    doPrint = True
    if doPrint and np.mod((k + 1 - k0MPC) * DT, 3600 / 10.0) == 0:
        print "====================================================="
        print 'Simulation time:', k * DT / 3600.0, 'hours'
        print 'PredictedMeasurmentError', (y - yP)  # /KF.measurementScaling
        print 'OptimControl', uOpt  # /scaleU
        print 'Objective', ca.DMatrix([Obj_k, Obj_total])
        if ULog is not None:
            print 'control DIFF', ULog['ct']
            print 'control USAT', ULog['uNS']
            print 'control cost', ULog['cost']

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
            Objh_NMPC.append(Obj_k)
            ObjTotalh.append(Obj_total)
        uh.append(uk)
        yPh.append(yP)
        xh.append(x_hat)
        yh.append(y)
        Juh_NMPC.append(Ju)
        simTime.append(k * DT)

toc = time.time()
computation_time = toc - tic
print 'computation time: ', computation_time

execfile('SavedResults\\plotCurves_model.py')