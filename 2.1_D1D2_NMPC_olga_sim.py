# -*- coding: utf-8 -*-
import sys
import os
import casadi as Ca
import control
import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as NP
from scipy.interpolate import interp1d
import ConfigParser

_path = os.path.dirname(os.path.realpath(__file__))

os.chdir(_path)
sys.path.append(_path+'\\NMPC')
sys.path.append(_path+'\\modelCompiler')
sys.path.append(_path+'\\StateEstimator')
sys.path.append(_path+'\\DaeModel')
sys.path.append(_path+'\\Olga_sim')

import modelCompiler
import SingleShooting
import MultipleShooting
import StateEstimator
import DaeModel
import Olga_sim as Olga

model = 'OlgaNet.ServerDemo.'
# serverIP = 'olga.itk.ntnu.no'
serverIP = 'localhost'
NP.set_printoptions(precision=3)
reload(StateEstimator)
reload(SingleShooting)
reload(MultipleShooting)
reload(Olga)
reload(DaeModel)

measurementTagsModelica = ['net.w1.w_L_out', 'net.w1.w_G_out', 'net.w2.w_L_out', 'net.w2.w_G_out', 'net.w1.P_r_t',
                           'net.w2.P_r_t', 'net.w1.deltaP_valve', 'net.w2.deltaP_valve', 'net.p.deltaP_valve']
measurementTagsOlga = ['WHD1.GLTHL', 'WHD1.GG', 'WHD2.GLTHL', 'WHD2.GG', 'WHD1.PT', 'WHD2.PT', 'Prod_Choke_1.VALVDP',
                       'Prod_Choke_2.VALVDP', 'Valve.VALVDP']

measurementsMonitorModelica = ['net.w1.fin.p', 'net.w2.fin.p', 'net.w1.P_bh', 'net.w2.P_bh', 'net.p.fin.p',
                               'net.p.w_L_out', 'net.p.w_G_out', 'net.PCA1.u', 'net.PCA2.u', 'net.PCINL.u',
                               'Obj', 'net.w1.P_res', 'net.w2.P_res']

measurementsMonitorOlga = ['ANN1.PT', 'ANN2.PT', 'BTH1.PT', 'BTH2.PT', 'PC-INL.MEASVAR',
                           'TOP.GLTHL', 'TOP.GG', 'PC-A1.CONTR', 'PC-A2.CONTR', 'PC-INL.CONTR',
                            'Cost.CONTR', 'IPR_1.RESPRESSURE', 'IPR_2.RESPRESSURE']

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
ocp = modelCompiler.getOCP('ocpGaslift_D1D2', models_list)

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

x0, z0, y_0 = daeModel.findSteadyState(u, None, None, 0, consList=['net.w1.z1', 'net.w2.z1', 'net.p.z'])

uOpt, xOpt, zOpt, yOpt = daeModel.findOptimalPoint(u, x0=x0, z0=z0, simCount=0,
                                                   consList=['net.w1.z1', 'net.w2.z1', 'net.p.z'])

x0 = xOpt
u0 = uOpt
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
R0 = NP.diag([1, 1, 1, 1, 1, 1, 1, 1, 1]) * 100
#R0 = NP.diag([100, 100, 100, 100, 1e5, 1e5, 1e4, 1e4, 1e4])
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

Config = ConfigParser.ConfigParser()
Config.read('config.ini')

MPC = True

if MPC:
    #    controlCost = 10000*NP.diag([1,1,1,1,1])
    #    finalStateCost = 0
    #    ss = SingleShooting.SingleShooting()  ## the steady state optimizer is implemented here.
    #    ss.defineOCP(ocp,DT=DTMPC,controlCost=controlCost,xOpt=xOpt,uOpt=uOpt,finalStateCost=finalStateCost)
    if Config.getboolean('MultipleShooting', 'MS'):
        solver = MultipleShooting.MultipleShooting()
    else:
        solver = SingleShooting.SingleShooting()

    deltaUCons = NP.array([0.02, 0.02, 1.0e3, 2.0e3, 2.0e3]) / DTMPC  ## given in SI

    ## check consistency
    if deltaUCons.size != 0:
        assert deltaUCons.size == ocp.u.size()
    controlCost = NP.diag([10.0, 10.0, 1.0, 1.0, 1.0])
    finalStateCost = 0.0
    solver.defineOCP(ocp, DT=DTMPC, controlCost=controlCost, xOpt=x0, uOpt=u0, finalStateCost=finalStateCost,
                     deltaUCons=deltaUCons)

# flush print buffer
print '----'
sys.stdout.flush()

## Start Simulator
openLoop = False
controlOlga = True
NIT = int(NP.ceil((3600 * 55 + 100)/DT))
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
simTime = []
Juh_NMPC = []
monitoredOlgah = []
monitoredModelicah = []
Objh_NMPC = []
ObjTotalh = []
xh.append(Ca.DMatrix(x0))
ULog = None

monitor = True
if monitor:
    measurementsEst = Ca.vertcat([ocp.variable(varName).v for varName in measurementsMonitorModelica])
    measure_funcEst = ocp.beq(measurementsEst)
    HEst = Ca.SXFunction(Ca.daeIn(x=ocp.x, z=ocp.z, p=ocp.u, t=ocp.t), [measure_funcEst])
    HEst.init()

extendedKalman = True
# u0 =  Ca.vertcat([ocp.variable(xit.getName()).initialGuess.getValue() for xit in ocp.u])
# x0,z0,y0 = daeModel.findSteadyState(u0,None,None,0,consList=['net.w1.z1','net.w2.z1','net.p.z'])


KF.setState(x0, z0)
x_hat = KF.getState()
z_hat = Ca.vertcat([ocp.variable(ocp.z[k].getName()).start for k in range(ocp.z.size())])

uk = u0

nX = ocp.x.size()
xMin = Ca.DMatrix(NP.zeros(nX))
xMax = Ca.DMatrix(NP.zeros(nX))
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
    print 'Obj. Value:', objValue/(3600*16)
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
        A_inv = NP.linalg.inv(-A)
        Ju = NP.array(C * (A_inv * B) + D)[0]

if openLoop or controlOlga:
    if controlOlga and not openLoop:
        y0Olga = Olga.OLGA_read_tags(da, measurementTagsOlga)

P_res1 = 160.0
P_res2 = 170.0
Obj_total = 0.0
dP_res = 5.0 / (10 * 3600)
gasPrice = ocp.variable('gasPrice').start

k0MPC = 5*3600 / DT
for k in range(0, NIT+1):
    sys.stdout.flush()

    if (k > 3600 * 10 / DT) and (k <= 3600 * 20 / DT):
        KF.R = R0 - 0.0275 * (k - 10*3600 / DT) * NP.diag(NP.ones(len(measurementTagsOlga)))

    if (k > 3600 * 10 / DT) and (k <= 3600 * 20 / DT):
        P_res1 -= dP_res * DT
    if (k > 3600 * 30 / DT) and (k <= 3600 * 40 / DT):
        P_res2 -= dP_res * DT

    if (NP.mod((k - k0MPC) * DT, DTMPC) == 0) and (k >= k0MPC):
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
        print 'Obj. Value:', objValue/(3600*16)
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

            A_inv = NP.linalg.inv(-A)
            Ju = NP.array(C * (A_inv * B) + D)[0]
        Olga.OLGA_write_tags(da, dict(zip(controlTagsOlga, NP.array(uk))))

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

        yMonitor[-2] = P_res1
        yMonitor[-1] = P_res2

        Obj_k = gasPrice * (uk[0] + uk[1]) - (y[0] + y[2])
        Obj_total = Obj_total + Obj_k

        doPrint = True
        if doPrint and NP.mod((k + 1 - k0MPC) * DT, DTMPC / 10.0) == 0:
            print "====================================================="
            print 'Simulation time:', (k + 1) * DT / (3600.0), 'hours'
            print 'PredictedMeasurmentError', (y - yP)  # /KF.measurementScaling
            print 'OptimControl', uOpt  # /scaleU
            print 'Objective', Ca.DMatrix([Obj_k, Obj_total])
            # print 'steadyStates',xOpt/scaleX
            #print '      states', x_hat
            # print 'steadyAlgStates',zOpt/scaleZ
            # print 'Alg      States',z_hat/scaleZ
            print 'Monitor    Olga', Ca.DMatrix(yMonitor)  # /monitorScale
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
        simTime.append(k*DT)

execfile('plotCurves_struc2.py')
execfile('SaveSimData_NMPC_olga.py')
