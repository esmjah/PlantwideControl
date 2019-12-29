# -*- coding: utf-8 -*-

# Authors: Andres Codas, Esmaeil Jahanshahi
# Norwegian University of Science and Technology

# This simulation requires the Olga Simulator. The Olga OPC Server must be running as Administrator user.
# The Olga case for this simulation is provider in the repository (\olgaModels\NetworkOPC\Network_Nominal.opi)
# This script runs NMPC in when process is at nominal conditions (no disturbances)
# Simulations starts from an initial point where the gas injection rate is 1 kg/s and all valves are 50% open
# Optimization steers the process to the steady-state optimal point
# The Modelica model "ocpGaslift_Nominal.mop" defines the Optimal Contrtrol Problem (OPC)

import sys
import os
import casadi as ca
import scipy.io as sio
import control
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import ConfigParser
import time

tic = time.time()

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
                               'net.w1.GOR_m_t', 'net.w2.GOR_m_t']
measurementsMonitorOlga = ['ANN1.PT', 'ANN2.PT', 'BTH1.PT', 'BTH2.PT', 'TOP.GLTHL', 'TOP.GG', 'PC-A1.CONTR',
                           'PC-A2.CONTR', 'PC-INL.CONTR', 'PC-INL.MEASVAR', 'GOR1.CONTR', 'GOR2.CONTR']

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
DTMPC = 600  ## for the MPC algorithm  Please always choose one to be multiple of the other
prediction_horizon = 10 * 3600

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

Config = ConfigParser.ConfigParser()
Config.read('config.ini')

MPC = True

if MPC:
    #  controlCost = 10000*np.diag([1,1,1,1,1])
    #  finalStateCost = 0
    #  ss = SingleShooting.SingleShooting()  ## the steady state optimizer is implemented here.
    #  ss.defineOCP(ocp,DT=DTMPC,controlCost=controlCost,xOpt=xOpt,uOpt=uOpt,finalStateCost=finalStateCost)
    if Config.getboolean('MultipleShooting', 'MS'):
        solver = MultipleShooting.MultipleShooting()
    else:
        solver = SingleShooting.SingleShooting()
    deltaUCons = np.array([0.015, 0.015, 1.0e3, 5.0e3, 5.0e3]) / DTMPC  # given in SI

    # check consistency
    if deltaUCons.size != 0:
        assert deltaUCons.size == ocp.u.size()

    controlCost = 2e3 * np.diag([10.0, 10.0, 1.0, 1.0, 1.0])
    finalStateCost = 0.0
    solver.defineOCP(ocp, DT=DTMPC, controlCost=controlCost, xOpt=x0, uOpt=u, finalStateCost=finalStateCost,
                     deltaUCons=deltaUCons)

# flush print buffer
print '----'
sys.stdout.flush()

# Start Simulator
openLoop = False
controlOlga = True
NIT = int(np.ceil(3600 * 36 / DT))
xSh = []
if not openLoop:
    if controlOlga:
        da = Olga.OLGA_connect(serverIP=serverIP, modelPrefix=model)
        # Olga.OLGA_restart_time(da,True)
    else:
        xSim = np.copy(x0) + ca.vertcat([np.random.normal() * scaleX[k] / 10 for k in range(ocp.x.size())])
        zSim = np.copy(z0)
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

        NIT = min(int(np.floor(yh[-1][OlgaTagIndex["TIME"]] / DT)), NIT)

        OlgaF = {}
        for k in olgaTags:
            OlgaF[k] = interp1d([yi[OlgaTagIndex['TIME']] for yi in yh], [yi[OlgaTagIndex[k]] for yi in yh])

    else:
        mh = []
        time = []
        x0k = x0 + [np.random.normal() * scaleX[k] / 10 for k in range(ocp.x.size())]
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
Juh_NMPC = []
monitoredOlgah = []
monitoredModelicah = []
Objh_NMPC = []
simTime = []
ObjTotalh = []
xh.append(ca.DMatrix(x0))
ULog = None

monitor = True
if monitor:
    measurementsEst = ca.vertcat([ocp.variable(varName).v for varName in measurementsMonitorModelica])
    measure_funcEst = ocp.beq(measurementsEst)
    HEst = ca.SXFunction(ca.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t), [measure_funcEst])
    HEst.init()


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

if MPC:
    uk, objValue, stats = solver.solveProblem(x0=x_hat, u0=u)
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

KalmanStopk = 2 * NIT

Obj_total = 0.0
gasPrice = ocp.variable('gasPrice').start

uk = u
k0 = int(0)
k0MPC = 3600 / DT
for k in range(k0, NIT):
    sys.stdout.flush()

    if (k > 3600 * 1 / DT) and (k <= 3600 * 11 / DT):
        KF.R = R0 - 0.027 * (k - 3600 / DT) * np.diag(np.ones(len(measurementTagsOlga)))

    if (np.mod((k - k0MPC) * DT, DTMPC) == 0) and (k >= k0MPC):
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
        Olga.OLGA_write_tags(da, dict(zip(controlTagsOlga, np.array(uk))))

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

        Obj_k = gasPrice * (uk[0] + uk[1]) - (y[0] + y[2])
        Obj_total = Obj_total + Obj_k

        doPrint = True
        if doPrint and np.mod((k + 1 - k0MPC) * DT, DTMPC / 10.0) == 0:
            print "====================================================="
            print 'Simulation time:', (k + 1) * DT / (3600.0), 'hours'
            print 'PredictedMeasurmentError', (y - yP)  # /KF.measurementScaling
            print 'Control input', uk  # /scaleU
            print 'Objective', ca.DMatrix([Obj_k, Obj_total])
            # print 'steadyStates',xOpt/scaleX
            print '      states', x_hat
            # print 'steadyAlgStates',zOpt/scaleZ
            # print 'Alg      States',z_hat/scaleZ
            print 'Monitor    Olga', ca.DMatrix(yMonitor)  # /monitorScale
            print 'MonitorModelica', HEst.getOutput()  # /monitorScale
            sys.stdout.flush()

    if k % 10 == 0:
        if monitor:
            monitoredOlgah.append(yMonitor)
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

execfile('SavedResults\\plotCurves_olga.py')
execfile('SavedResults\\SaveSimData_NMPC_olga.py')
